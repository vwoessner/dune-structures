from dune.codegen.options import get_option
from dune.codegen.ufl.execution import *
from dune.structures.codegen import (
    UFLPhysicalParameter,
    UFLMaterialLawIndex,
    UFLPrestress,
)
from dune.testtools.parametertree.parser import parse_ini_file
from ufl.cell import simplex


def _elasticity_form_impl(u, v, cell, materials):
    law_index = UFLMaterialLawIndex(cell)

    def material_form(material):
        piola = material.first_piola(u)

        # Add active prestressed (isotropic) material
        stress = piola + UFLPrestress(cell)

        dxm = dx(subdomain_data=law_index, subdomain_id=material.id)
        return inner(stress, grad(v)) * dxm

    # The final form is a sum of all possible materials
    form = sum((material_form(m) for m in materials), Form([]))

    return form


def elasticity_form(materials, degree=1, dim=3):
    # Define cell
    cell = simplex(dim)

    # Setup finite elements
    element = VectorElement("CG", cell, degree)
    u = TrialFunction(element)
    v = TestFunction(element)

    form = _elasticity_form_impl(u, v, cell, materials)

    # Add some body force (source term)
    force = Coefficient(element, cargo={"name": "Force"})
    form = form - inner(force, v) * dx

    # Add some surface traction (Neumann boundary condition)
    traction = Coefficient(element, cargo={"name": "Traction"})
    form = form - inner(traction, v) * ds

    return form


def elastodynamics_form(materials, degree=1, dim=3):
    # Define cell
    cell = simplex(dim)

    # Setup finite elements
    delement = VectorElement("CG", cell, degree)
    element = MixedElement(delement, delement)

    u0, u1 = TrialFunctions(element)
    v0, v1 = TestFunctions(element)

    mass = (inner(u0, v1) + inner(u1, v0)) * dx
    r = _elasticity_form_impl(u0, v0, cell, materials) - inner(u1, v1) * dx

    return r, mass


def quasistatic_mass_form(degree=1, dim=3):
    # Define cell
    cell = simplex(dim)

    # Setup finite elements
    element = VectorElement("CG", cell, degree)
    u = TrialFunction(element)
    v = TestFunction(element)

    # Later attach some meaning to this value
    attenuation = 1e-6

    return attenuation * inner(u, v) * dx