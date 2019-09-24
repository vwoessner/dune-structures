""" An abstraction for material laws for nonlinear mechanics.
This is inspired by a retired Simula project I came along:
https://bazaar.launchpad.net/~cbc-core/cbc.solve/main/files/head:/cbc/twist

New material models are defined in terms of a strain energy functional.
From that, the weak formulation is automatically derived.
"""

from dune.codegen.ufl.execution import *
from dune.structures.codegen import UFLPhysicalParameter


class MaterialLawBase(object):
    #
    # The interface methods that need to be implemented for each new material law
    #
    def strain_energy(self, u):
        raise NotImplementedError
    
    @property
    def strain_measure_type(self):
        raise NotImplementedError

    @property
    def param_names(self):
        raise NotImplementedError

    #
    # Mapping physical parameters to UFL objects. Right now this defaults to an
    # implementation that supports runtime heterogeneity.
    #
    def ufl_parameters(self, u):
        cell = u.ufl_element().cell()
        return [UFLPhysicalParameter(p, cell) for p in self.param_names]

    def strain_measure(self, u):
        if not hasattr(self, self.strain_measure_type):
            raise NotImplementedError("Strain measure '{}' not known".format(self.strain_measure_type))
        return getattr(self, self.strain_measure_type)(u)

    #
    # Piola-Kirchhoff tensors as UFL expressions
    #
    def first_piola(self, u):
        S = self.second_piola(u)

        if self.strain_measure_type == "infinitesimal_strain":
            return S
        else:
            F = Identity(3) + grad(u)
            return dot(F, S)

    def second_piola(self, u):
        energy = self.strain_energy(u)
        strain = self.strain_measure(u)

        if self.strain_measure_type == "infinitesimal_strain":
            return diff(energy, strain)
        else:
            raise NotImplementedError

    #
    # Implementations of strain measures to be used by energy functionals
    #
    def infinitesimal_strain(self, u):
        from ufl.classes import Variable, Label
        return Variable(0.5 * (grad(u) + grad(u).T), label=Label(42))


class LinearMaterial(MaterialLawBase):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def strain_measure_type(self):
        return "infinitesimal_strain"

    def strain_energy(self, u):
        epsilon = self.strain_measure(u)
        mu, lmbda = self.ufl_parameters(u)
        return 0.5 * lmbda * (tr(epsilon)**2) + mu*tr(epsilon * epsilon)
