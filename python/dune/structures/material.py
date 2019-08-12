""" A material class that wraps our C++ material concept """

from dune.codegen.generation import (class_member,
                                     constructor_parameter,
                                     initializer_list,
                                     function_mangler,
                                     template_parameter,
                                     globalarg,
                                     )
from dune.codegen.pdelab.geometry import (enforce_boundary_restriction,
                                          name_cell,
                                          world_dimension,
                                          )
from dune.codegen.loopy.target import dtype_floatingpoint
from dune.codegen.tools import maybe_wrap_subscript
import pymbolic.primitives as prim
import loopy as lp
import numpy as np
import ufl


class ElasticMaterialBase(object):
    def first_lame(self):
        raise NotImplementedError
    
    def second_lame(self):
        raise NotImplementedError


class HomogeneousMaterial(ElasticMaterialBase):
    def __init__(self, lame1, lame2):
        self.lame1 = lame1
        self.lame2 = lame2

    def first_lame(self):
        return self.lame1
    
    def second_lame(self):
        return self.lame2


class HeterogeneousMaterial(ElasticMaterialBase):
    def first_lame(self):
        return UFLPhysicalParameter("first_lame")

    def second_lame(self):
        return UFLPhysicalParameter("second_lame")


class UFLPhysicalParameter(ufl.coefficient.Coefficient):
    def __init__(self, param):
        self.param = param
#         FE = ufl.FiniteElement("Real", 'tetrahedron', 0)
        FE = ufl.FiniteElement("DG", 'tetrahedron', 0)
        ufl.coefficient.Coefficient.__init__(self, FE)
    
    def _ufl_expr_reconstruct_(self):
        return UFLPhysicalParameter(self.param)

    def visit(self, visitor):
        restriction = enforce_boundary_restriction(visitor)
        cell = name_cell(restriction)
        globalarg(cell, dtype=lp.types.to_loopy_type(np.dtype(str)))
        cell = prim.Variable(cell)
        # The following is a terrible hack until we get the function interface
        # in loopy: We need to pass the quadrature points coordinate by coordinate
        # and piece the information back together on the receiving side.
        qp = tuple(maybe_wrap_subscript(visitor.to_cell(visitor.quadrature_position()), i) for i in range(world_dimension()))
        return prim.Call(LoopyPhysicalParameter(self.param), (cell,) + qp)


class LoopyPhysicalParameter(lp.symbolic.FunctionIdentifier):
    def __init__(self, param):
        self.param = param

    def __getinitargs__(self):
        return (self.param,)
    
    @property
    def name(self):
        material_class = name_material_class()
        return "{}.{}".format(material_class, self.param) 


@function_mangler
def lame_parameter_function_mangler(knl, func, arg_dtypes):
    if isinstance(func, LoopyPhysicalParameter):
        return lp.CallMangleInfo(func.name,
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),),
                                 (lp.types.to_loopy_type(np.dtype(str)),) +
                                 (lp.types.to_loopy_type(dtype_floatingpoint()),) * world_dimension()
                                )

def name_material_class():
    name = "material"
    define_material_class(name)
    return name


@class_member(classtag="operator")
def define_material_class(name):
    tp = name_material_template_parameter()
    constructor_parameter("const {}&".format(tp), name, classtag="operator")
    initializer_list(name, [name], classtag="operator")
    return "{} {};".format(tp, name)


@template_parameter(classtag="operator")
def name_material_template_parameter():
    return "MATERIAL"
