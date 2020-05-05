""" Code generator for some differential geometry quantities.

This is not embedded into the dune-codegen code generation pipeline.
Instead, it is meant to be run by hand and the resulting expressions
are copied and pasted into the 1D fibre operator.

Currently, this script needs two adjustments in codegen to be run:
* In python/ufl/transformations/__init__.py  the transformation printing needs to be disabled
* In python/ufl/modified_terminals.py the assertion in the grad handler needs to be disabled
I am too lazy right now to find proper upstream fixes for these minor issues.
"""

from dune.codegen.ufl.execution import *
from dune.codegen.pdelab.restriction import restricted_name
from dune.codegen.pdelab.localoperator import GenericAccumulationMixin
from ufl.algorithms import expand_derivatives, extract_arguments
from ufl.algorithms.apply_restrictions import apply_restrictions
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

    # Tangential and normal vector
    t = Coefficient(FE, cargo={"name": "t", "diff": 0, "restriction": 0})
    n = perp(t)

    cell_quantities = {
        "dtut" : inner(t, grad(inner(u, t))),
        "dtvt" : inner(t, grad(inner(v, t))),
        "dt2un" : inner(t, grad(inner(t, grad(inner(u, n))))),
        "dt2vn" : inner(t, grad(inner(t, grad(inner(v, n))))),
    }

    facet_quantities = {
        "sk_dt2un" : avg(inner(t, grad(inner(t, grad(inner(u, n)))))),
        "dt2vn_n" : inner(t, grad(inner(t, grad(inner(v('+'), n))))),
        "dt2vn_s" : inner(t, grad(inner(t, grad(inner(v('-'), n))))),
        "sk_dtun" : jump(inner(t, grad(inner(u, n)))),
        "dtvn_n" : inner(t, grad(inner(v('+'), n))),
        "dtvn_s" : inner(t, grad(inner(v('-'), n))),
    }

    def print_code(name, expr, measure):
        # Do the preprocessing
        expr = apply_function_pullbacks(expr)
        expr = expand_derivatives(expr)
        expr = apply_algebra_lowering(expr)
        expr = remove_complex_nodes(expr)

        # This weird check would go away if we reproduced more of the UFL preprocessing here...
        if name.startswith('sk'):
            expr = apply_restrictions(expr)

        expr = pushdown_indexed(expr)

        print("\n\n\nQuantity: {}".format(name))
        for snippet in ufl_to_code(expr, measure):
            print(snippet + '\n\n')

    for name, expr in cell_quantities.items():
        print_code(name, expr, "cell")

    for name, expr in facet_quantities.items():
        print_code(name, expr, "interior_facet")


class AdHocVisitor(GenericAccumulationMixin, UFL2LoopyVisitor):
    def __init__(self, measure):
        UFL2LoopyVisitor.__init__(self, "cell", {})
        self.grad_count = 0
        self.measure = measure

    def __call__(self, o, do_split=False):
        if len(extract_arguments(o)) == 0:
            return [self._call(o, False)]
        else:
            rets = []
            for info in self.list_accumulation_infos(o):
                self.current_info = info
                expr = self._call(o, False)
                if expr != 0:
                    rets.append(expr)
            return rets

    def lfs_inames(self, *args):
        return []

    def argument(self, o):
        info = self.get_accumulation_info(o)
        if info != self.current_info[0]:
            self.indices = None
            return 0

        name = restricted_name("basis", self.restriction)

        if self.grad_count == 0:
            name = "{}.function".format(name)
        elif self.grad_count == 1:
            name = "{}.jacobian".format(name)
        elif self.grad_count == 2:
            name = "{}.hessian".format(name)
        else:
            raise NotImplementedError

        indices = self.indices[1:]
        self.indices = None
        return prim.Call(prim.Variable(name), (prim.Variable("i"),) + indices)

    def coefficient(self, o):
        name = restricted_name(o.cargo["name"], o.cargo.get("restriction", self.restriction))

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

        return prim.Subscript(prim.Variable(restricted_name("jit", self.restriction)), (j, i))


class IndexNester(IdentityMapper):
    def map_subscript(self, expr):
        if len(expr.index) > 1:
            return self.rec(prim.Subscript(prim.Subscript(expr.aggregate, expr.index[:1]), expr.index[1:]))
        else:
            return expr

def ufl_to_code(expr, measure):
    from dune.codegen.generation import global_context
    with global_context(integral_type="cell", form_identifier="foo"):
        visitor = AdHocVisitor(measure)
        nester = IndexNester()
        from pymbolic.mapper.c_code import CCodeMapper
        ccm = CCodeMapper()
        exprs = visitor(expr)
        exprs = [nester(e) for e in exprs]
        return [ccm(e) for e in exprs]
