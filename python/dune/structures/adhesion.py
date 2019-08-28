""" Models for focal adhesion and their UFL implementation """
from dune.codegen.ufl.execution import *
from pytools import ImmutableRecord


class DeshpandeAdhesionParameters(ImmutableRecord):
    def __init__(self):
        ImmutableRecord.__init__(self,
            bond_stiffness=1.5e-5,
            peak_bond_length=1.3e-8,
            total_integrin_concentration=5e-9,
            internal_energy_difference=2.14e-20,
            temperature=300.0,
            boltzmann_constant=1.38e-23,
        )

def deshpande_adhesion(form, where):
    # Extract test and trial function
    displacement = form.coefficients()[0]
    test = form.arguments()[0]

    # Instantiate the parameter object
    params = DeshpandeAdhesionParameters()

    # Calculate the potential
    peak = params.peak_bond_length
    length = sqrt(displacement[0] ** 2 + displacement[1] ** 2)
    potential = conditional(length < peak,
                            length ** 2,
                            conditional(length < 2.0 * peak,
                                        -peak ** 2 + 2.0 * length * peak - length ** 2,
                                        peak ** 2
                                        )
                            )
    potential = params.bond_stiffness * potential

    # Calculate the force - derivative of the potential
    force = tuple(
        conditional(length < peak,
                    2.0 * displacement[i],
                    conditional(length < 2.0 * peak,
                                2.0 * peak * displacement[i] / length - displacement[i],
                                0.0)
                    )
        for i in range(2)
        )
    force = tuple(params.bond_stiffness * f for f in force)

    # Calculate high affinity integrin concentration
    sumforce = force[0] * displacement[0] + force[1] * displacement[1]
    high_concentration = params.total_integrin_concentration / (exp(-(params.internal_energy_difference + potential - sumforce) / (params.temperature * params.boltzmann_constant)) + 1)

    # Finally calculate the traction
    traction = tuple(-high_concentration * force[i] for i in range(2)) + (0.0,)

    # Alter the form
    global ds
    ds = ds(subdomain_id=42, subdomain_data=conditional(where, 42, 0))
    return form - inner(as_vector(traction), test) * ds
