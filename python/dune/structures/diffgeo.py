""" Code generator for some differential geometry quantities.

This is not embedded into the dune-codegen code generation pipeline.
Instead, it is meant to be run by hand and the resulting expressions
are copied and pasted into the 1D fibre operator.
"""

from dune.codegen.ufl.execution import *
from ufl.algorithms import expand_derivatives
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.algorithms.remove_complex_nodes import remove_complex_nodes
from ufl.algorithms.apply_function_pullbacks import apply_function_pullbacks
from dune.codegen.ufl.transformations.indexpushdown import pushdown_indexed
from dune.codegen.ufl.visitor import UFL2LoopyVisitor
from pymbolic.mapper import IdentityMapper
import pymbolic.primitives as prim



def generate_tangential_derivatives():
    # The 2 here is the polynomial order of the finite element
    FE = VectorElement("CG", triangle, 2)
    u = Coefficient(FE, cargo={"name": "u"})
    v = TestFunction(FE)
    t = Coefficient(FE, cargo={"name": "t", "diff": 0})
    n = Coefficient(FE, cargo={"name": "n", "diff": 0})

    quantities = {
        "dtut" : inner(t, grad(inner(u, t))),
        "dtvt" : inner(t, grad(inner(v, t))),
        "dt2un" : inner(t, grad(inner(t, grad(inner(u, n))))),
        "dt2vn" : inner(t, grad(inner(t, grad(inner(v, n))))),
    }

    for name, expr in quantities.items():
        # Do the preprocessing
        expr = apply_function_pullbacks(expr)
        expr = expand_derivatives(expr)
        expr = apply_algebra_lowering(expr)
        expr = remove_complex_nodes(expr)
        expr = pushdown_indexed(expr)
    
        print("\n\nThe expression for the quantity {} is:".format(name))
        print(ufl_to_code(expr))


class AdHocVisitor(UFL2LoopyVisitor):
    def __init__(self):
        UFL2LoopyVisitor.__init__(self, "cell", {})
        self.grad_count = 0
        
    def argument(self, o):
        name = "phi"
        if self.grad_count > 0:
            name = "d{}{}".format(self.grad_count, name)
        return prim.Subscript(prim.Variable(name), (prim.Variable("i"), 0))

    def coefficient(self, o):
        name = o.cargo["name"]
        if self.grad_count > 0:
            if "diff" in o.cargo:
                self.indices = None
                return o.cargo["diff"]
            name = "d{}{}".format(self.grad_count, name)
        return prim.Variable(name)
    
    def reference_grad(self, o):
        self.grad_count = self.grad_count + 1
        ret = self.call(o.ufl_operands[0])
        self.grad_count = self.grad_count - 1
        return ret
    
    def jacobian_inverse(self, o):
        i, j = self.indices
        self.indices = None
        
        return prim.Subscript(prim.Variable("jit"), (j, i))


class IndexNester(IdentityMapper):
    def map_subscript(self, expr):
        if len(expr.index) > 1:
            return self.rec(prim.Subscript(prim.Subscript(expr.aggregate, expr.index[:1]), expr.index[1:]))
        else:
            return expr

def ufl_to_code(expr):
    from dune.codegen.generation import global_context
    with global_context(integral_type="cell", form_identifier="foo"):
        visitor = AdHocVisitor()
        nester = IndexNester()
        from pymbolic.mapper.c_code import CCodeMapper
        ccm = CCodeMapper()
        
        expr = visitor(expr)
        
        # Not worth it right now, maybe for some other expressions
        from dune.codegen.sympy import simplify_pymbolic_expression
        expr = simplify_pymbolic_expression(expr)
         
        expr = nester(expr)
        return ccm(expr)