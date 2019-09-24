from dune.codegen.ufl.execution import *
from dune.structures.codegen import UFLPhysicalParameter


def elasticity_form(material):
    # Define cell
    cell = tetrahedron

    # Setup finite elements
    element = VectorElement("CG", cell, 1)
    u = TrialFunction(element)
    v = TestFunction(element)

    piola = material.first_piola(u)

    # Add active prestressed (isotropic) material
    stress = piola + UFLPhysicalParameter("pretension", cell) * Identity(3)

    return inner(stress, grad(v)) * dx
