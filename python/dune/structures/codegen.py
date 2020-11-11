""" Extensions for the dune-codegen code generator.
These do the interfacing to our C++ abstraction of heterogeneous materials.
"""

from dune.codegen.ufl.execution import Coefficient
from dune.codegen.ufl.modified_terminals import Restriction
from dune.codegen.generation import (base_class,
                                     class_member,
                                     constructor_parameter,
                                     delete_cache_items,
                                     hook,
                                     include_file,
                                     initializer_list,
                                     instruction,
                                     function_mangler,
                                     template_parameter,
                                     temporary_variable,
                                     )
from dune.codegen.pdelab.geometry import (enforce_boundary_restriction,
                                          name_element_geometry_wrapper,
                                          world_dimension,
                                          )
from dune.codegen.pdelab.localoperator import lop_template_ansatz_gfs, lop_template_test_gfs
from dune.codegen.loopy.target import dtype_floatingpoint, type_floatingpoint
from dune.codegen.tools import maybe_wrap_subscript
import dune.codegen
import pymbolic.primitives as prim
import loopy as lp
import numpy as np
import ufl


class UFLPhysicalParameter(Coefficient):
    def __init__(self, index, cell):
        self.index = index
        self.cell = cell
        FE = ufl.FiniteElement("DG", cell, 0)
        Coefficient.__init__(self, FE)

    def _ufl_expr_reconstruct_(self):
        return UFLPhysicalParameter(self.index, self.cell)

    def visit(self, visitor):
        restriction = enforce_boundary_restriction(visitor)
        cell = prim.Call(LoopyEntity(restriction), ())
        # The following is a terrible hack until we get the function interface
        # in loopy: We need to pass the quadrature points coordinate by coordinate
        # and piece the information back together on the receiving side.
        qp = tuple(maybe_wrap_subscript(visitor.to_cell(visitor.quadrature_position()), i) for i in range(world_dimension()))
        return prim.Call(LoopyPhysicalParameter(), (cell, self.index) + qp)


class UFLMaterialLawIndex(Coefficient):
    def __init__(self, cell):
        self.cell = cell
        FE = ufl.FiniteElement("DG", cell, 0)
        Coefficient.__init__(self, FE)

    def _ufl_expr_reconstruct_(self):
        return UFLMaterialLawIndex(self.cell)

    def visit(self, visitor):
        restriction = enforce_boundary_restriction(visitor)
        cell = prim.Call(LoopyEntity(restriction), ())
        # The following is a terrible hack until we get the function interface
        # in loopy: We need to pass the quadrature points coordinate by coordinate
        # and piece the information back together on the receiving side.
        return prim.Call(LoopyMaterialLawIndex(), (cell,))


class UFLPrestress(Coefficient):
    def __init__(self, cell):
        self.cell = cell
        dim = cell.topological_dimension()
        FE = ufl.TensorElement("DG", cell, 0, (dim, dim))
        Coefficient.__init__(self, FE)

    def _ufl_expr_reconstruct(self):
        return UFLPrestress(self.cell)

    def visit(self, visitor):
        restriction = enforce_boundary_restriction(visitor)

        # Get the entity as a string instead of a proper loopy call, because
        # we need to construct a C-string.
        from dune.codegen.pdelab.geometry import name_cell
        cell = name_cell(restriction)

        name = "prestress_eval"
        dim = self.cell.topological_dimension()
        temporary_variable(name,
                           shape=(dim, dim),
                           managed=False,
                           shape_impl=("fm",))
        material_class = name_material_class()

        instruction(code="{}->prestress({}, {}, {});".format(material_class, str(cell), str(visitor.quadrature_position()), name),
                    assignees=frozenset({name}),
                    within_inames=frozenset(visitor.quadrature_inames()),
                    )
        return prim.Variable(name)


class LoopyEntity(lp.symbolic.FunctionIdentifier):
    def __init__(self, restriction):
        # Too lazy to implement this for skeletons, given that we do not need it
        assert restriction == Restriction.NONE
        self.restriction = restriction

    def __getinitargs__(self):
        return (self.restriction,)

    @property
    def name(self):
        eg = name_element_geometry_wrapper()
        return "{}.entity".format(eg)


class LoopyPhysicalParameter(lp.symbolic.FunctionIdentifier):
    def __getinitargs__(self):
        return ()

    @property
    def name(self):
        material_class = name_material_class()
        return "{}->parameter_unrolled".format(material_class)


class LoopyMaterialLawIndex(lp.symbolic.FunctionIdentifier):
    def __getinitargs__(self):
        return ()

    @property
    def name(self):
        material_class = name_material_class()
        return "{}->material_law_index".format(material_class)


@function_mangler
def dune_structures_function_mangler(knl, func, arg_dtypes):
    if isinstance(func, LoopyPhysicalParameter):
        return lp.CallMangleInfo(func.name,
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),),
                                 (lp.types.to_loopy_type(np.dtype(str)),
                                  lp.types.to_loopy_type(np.int32)) +
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),) * world_dimension()
                                )

    if isinstance(func, LoopyMaterialLawIndex):
        return lp.CallMangleInfo(func.name,
                                 (lp.types.to_loopy_type(np.int32),),
                                 (lp.types.to_loopy_type(np.dtype(str)),)
                                 )

    if isinstance(func, LoopyEntity):
        return lp.CallMangleInfo(func.name, (lp.types.to_loopy_type(np.dtype(str)),), ())


def name_material_class():
    name = "material"
    define_material_class(name)
    return name


@class_member(classtag="operator")
def define_material_class(name):
    include_file("dune/structures/material.hh", filetag="operatorfile")
    _type = "std::shared_ptr<ElasticMaterialBase<typename {}::Traits::EntitySet, {}>>".format(lop_template_ansatz_gfs(), type_floatingpoint())
    constructor_parameter(_type, name, classtag="operator")
    initializer_list(name, [name], classtag="operator")
    setter_material_class(name)
    return "{} {};".format(_type, name)


@class_member(classtag="operator")
def setter_material_class(name):
    _type = "std::shared_ptr<ElasticMaterialBase<typename {}::Traits::EntitySet, {}>>".format(lop_template_ansatz_gfs(), type_floatingpoint())
    return ["void setMaterial({} {}_)".format(_type, name),
            "{",
            "  {} = {}_;".format(name, name),
            "}"
            ]


# Monkey-patch construct_signature to drop template parameters
def notemplate_construct_signature(types, args, name):
    for _type in types:
        using_from_baseclass(_type)
    func = "void {}({}) const override".format(name, ", ".join("{}{}& {}".format("const " if c else "", t, a) for t, (c, a) in zip(types, args)))
    return [func]


dune.codegen.pdelab.signatures.construct_signature = notemplate_construct_signature


@class_member(classtag="operator")
def using_from_baseclass(_type):
    gfsu = lop_template_ansatz_gfs()
    gfsv = lop_template_test_gfs()
    base = 'Dune::BlockLab::AbstractLocalOperatorInterface<{}, {}>'.format(gfsu, gfsv)
    return "using {} = typename {}::{};".format(_type, base, _type)


# Add the abstract base class for local operators from dune-structures
@hook('operator_class_extraction')
def baseclass_hook():
    gfsu = lop_template_ansatz_gfs()
    gfsv = lop_template_test_gfs()
    include_file('dune/blocklab/operators/virtualinterface.hh', filetag="operatorfile")
    delete_cache_items("baseclass")
    base_class('Dune::BlockLab::AbstractLocalOperatorInterface<{}, {}>'.format(gfsu, gfsv), classtag="operator")
