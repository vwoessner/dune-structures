# Make sure that monnkey patches are triggered
import dune.structures.codegen

from dune.structures.adhesion import deshpande_adhesion

from dune.structures.elasticity import (
    elasticity_form,
    elastodynamics_form,
    quasistatic_mass_form,
)

from dune.structures.material import (
    LinearMaterial,
    NeoHookeanMaterial,
    StVenantKirchhoffMaterial,
)
