from dune.codegen.options import get_option
from dune.codegen.ufl.execution import *
from dune.structures.codegen import UFLPhysicalParameter, UFLMaterialLawIndex
from dune.testtools.parametertree.parser import parse_ini_file


def elasticity_form(materials, force=None):
    # Apply defaults
    if force is None:
        force = as_vector([0.0, 0.0, 0.0])

    # Define cell
    cell = tetrahedron

    # Setup finite elements
    element = VectorElement("CG", cell, 1)
    u = TrialFunction(element)
    v = TestFunction(element)

    law_index = UFLMaterialLawIndex(cell)

    def material_form(material):
        piola = material.first_piola(u)
        # Add active prestressed (isotropic) material
        stress = piola + UFLPhysicalParameter("pretension", cell) * Identity(3)

        dxm = dx(subdomain_data=law_index, subdomain_id=material.id)
        return inner(stress, grad(v)) * dxm

    # The final form is a sum of all possible materials
    form = sum((material_form(m) for m in materials), Form([]))

    # This works around the stupid zero elimination bug
    if any(f != 0.0 for f in force):
        form = form - inner(force, v) * dx

    return form
