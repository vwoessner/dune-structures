from dune.codegen.ufl.execution import *
import ufl.classes


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

    stress = material.first_lame() * div(u) * Identity(3) + 2.0 * material.second_lame() * sym(grad(u))

    # Add pretension
    stress = stress + material.pretension() * Identity(3)

    # The linear elasticity form
    return (inner(stress, sym(grad(v))) - inner(force, v)) * dx
