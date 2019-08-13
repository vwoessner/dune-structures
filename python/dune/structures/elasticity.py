from dune.codegen.ufl.execution import *
import ufl.classes

class ElasticityBCType:
    NEUMANN = 0
    DIRICHLET = 1


def sym(expr):
    return 0.5 * (expr + expr.T)


def linear_elasticity_form(material,
                           force=None):
    # Apply defaults
    if force is None:
        force = as_vector([0.0, 0.0, 0.0])

    # Define cell and set it in the material class
    cell = tetrahedron
    material = material.with_cell(cell)

    # Setup finite elements
    element = VectorElement("CG", cell, 1)
    u = TrialFunction(element)
    v = TestFunction(element)

    # The linear elasticity form
    return (inner(material.first_lame() * div(u) * Identity(3) + 2.0 * material.second_lame() * sym(grad(u)), sym(grad(v))) - inner(force, v)) * dx
