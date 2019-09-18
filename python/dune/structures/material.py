""" A material class that wraps our C++ material concept """

from dune.codegen.ufl.modified_terminals import Restriction
from dune.codegen.generation import (class_member,
                                     constructor_parameter,
                                     include_file,
                                     initializer_list,
                                     function_mangler,
                                     template_parameter,
                                     )
from dune.codegen.pdelab.geometry import (enforce_boundary_restriction,
                                          name_element_geometry_wrapper,
                                          world_dimension,
                                          )
from dune.codegen.pdelab.localoperator import lop_template_ansatz_gfs
from dune.codegen.loopy.target import dtype_floatingpoint, type_floatingpoint
from dune.codegen.tools import maybe_wrap_subscript
import pymbolic.primitives as prim
import loopy as lp
import numpy as np
import ufl


class ElasticMaterialBase(object):
    def first_lame(self):
        """ Return the first lame coefficient (often lambda) """
        raise NotImplementedError
    
    def second_lame(self):
        """ Return the second lame coeffcient (often mu) """
        raise NotImplementedError

    def pretension(self):
        """ Return the pretension constant """
        raise NotImplementedError

    def with_cell(self, cell):
        """ Bind this material to a given cell - default is no op """
        return self


class HomogeneousMaterial(ElasticMaterialBase):
    def __init__(self, lame1, lame2, pretension=0.0):
        self.lame1 = lame1
        self.lame2 = lame2
        self.pretens = pretension

    def first_lame(self):
        return self.lame1

    def second_lame(self):
        return self.lame2

    def pretension(self):
        return self.pretens


class HeterogeneousMaterial(ElasticMaterialBase):
    def __init__(self, cell=None):
        self.cell = cell

    def with_cell(self, cell):
        return HeterogeneousMaterial(cell)

    def first_lame(self):
        return UFLPhysicalParameter("first_lame", self.cell)

    def second_lame(self):
        return UFLPhysicalParameter("second_lame", self.cell)

    def pretension(self):
        return UFLPhysicalParameter("pretension", self.cell)


class UFLPhysicalParameter(ufl.coefficient.Coefficient):
    def __init__(self, param, cell):
        self.param = param
        self.cell = cell
        FE = ufl.FiniteElement("DG", cell, 0)
        ufl.coefficient.Coefficient.__init__(self, FE)
    
    def _ufl_expr_reconstruct_(self):
        return UFLPhysicalParameter(self.param, self.cell)

    def visit(self, visitor):
        restriction = enforce_boundary_restriction(visitor)
        cell = prim.Call(LoopyEntity(restriction), ())
        # The following is a terrible hack until we get the function interface
        # in loopy: We need to pass the quadrature points coordinate by coordinate
        # and piece the information back together on the receiving side.
        qp = tuple(maybe_wrap_subscript(visitor.to_cell(visitor.quadrature_position()), i) for i in range(world_dimension()))
        return prim.Call(LoopyPhysicalParameter(self.param), (cell,) + qp)


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
    def __init__(self, param):
        self.param = param

    def __getinitargs__(self):
        return (self.param,)
    
    @property
    def name(self):
        material_class = name_material_class()
        return "{}->{}".format(material_class, self.param)


@function_mangler
def dune_structures_function_mangler(knl, func, arg_dtypes):
    if isinstance(func, LoopyPhysicalParameter):
        return lp.CallMangleInfo(func.name,
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),),
                                 (lp.types.to_loopy_type(np.dtype(str)),) +
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),) * world_dimension()
                                )

    if isinstance(func, LoopyEntity):
        return lp.CallMangleInfo(func.name, (lp.types.to_loopy_type(np.dtype(str)),), ())


def name_material_class():
    name = "material"
    define_material_class(name)
    return name


@class_member(classtag="operator")
def define_material_class(name):
    include_file("dune/structures/material.hh", filetag="operator")
    _type = "std::shared_ptr<ElasticMaterialBase<typename {}::Traits::EntitySet, {}>>".format(lop_template_ansatz_gfs(), type_floatingpoint())
    constructor_parameter(_type, name, classtag="operator")
    initializer_list(name, [name], classtag="operator")
    return "{} {};".format(_type, name)
