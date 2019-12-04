from dune.codegen.options import get_option
from dune.codegen.ufl.execution import *
from dune.structures.codegen import UFLPhysicalParameter, UFLMaterialLawIndex
from dune.testtools.parametertree.parser import parse_ini_file


def _elasticity_form_impl(u, v, cell, materials, force):
    law_index = UFLMaterialLawIndex(cell)

    def material_form(material):
        piola = material.first_piola(u)
        # Add active prestressed (isotropic) material
        stress = piola + UFLPhysicalParameter("pretension", cell) * Identity(3)

        dxm = dx(subdomain_data=law_index, subdomain_id=material.id)
        return inner(stress, grad(v)) * dxm

    # The final form is a sum of all possible materials
    form = sum((material_form(m) for m in materials), Form([]))

    # Maybe add some body force
    if force:
        f = Coefficient(element)
        form = form - inner(force, v) * dx

    return form



def elasticity_form(materials, force=False):
    # Define cell
    cell = tetrahedron

    # Setup finite elements
    element = VectorElement("CG", cell, 1)
    u = TrialFunction(element)
    v = TestFunction(element)

    return _elasticity_form_impl(u, v, cell, materials, force)


def elastodynamics_form(materials):
    # Define cell
    cell = tetrahedron

    # Setup finite elements
    delement = VectorElement("CG", cell, 1)
    element = MixedElement(delement, delement)

    u0, u1 = TrialFunctions(element)
    v0, v1 = TestFunctions(element)

    mass = (inner(u0, v1) + inner(u1, v0)) * dx
    r = _elasticity_form_impl(u0, v0, cell, materials, False) - inner(u1, v1) * dx

    return r, mass