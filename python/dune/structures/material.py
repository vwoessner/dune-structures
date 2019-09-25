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
    def param_names(self):
        raise NotImplementedError

    @property
    def id(self):
        raise NotImplementedError

    #
    # Mapping physical parameters to UFL objects. Right now this defaults to an
    # implementation that supports runtime heterogeneity.
    #
    def ufl_parameters(self, u):
        cell = u.ufl_element().cell()
        return [UFLPhysicalParameter(p, cell) for p in self.param_names]

    #
    # Some physical quantities for convenience and use in derived classes 
    #
    def first_piola(self, u):
        S = self.second_piola(u)
        F = self.deformation_gradient(u)
        return dot(F, S)

    def second_piola(self, u):
        raise NotImplementedError

    def deformation_gradient(self, u):
        return Identity(3) + grad(u)

    def infinitesimal_strain(self, u):
        from ufl.classes import Variable, Label
        return Variable(0.5 * (grad(u) + grad(u).T), label=Label(42))


class LinearMaterial(MaterialLawBase):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def id(self):
        return 0

    def deformation_gradient(self, u):
        return Identity(3)

    def strain_energy(self, u):
        epsilon = self.infinitesimal_strain(u)
        mu, lmbda = self.ufl_parameters(u)
        return 0.5 * lmbda * (tr(epsilon)**2) + mu*tr(epsilon * epsilon)

    def second_piola(self, u):
        return diff(self.strain_energy(u), self.infinitesimal_strain(u))


class NeoHookeanMaterial(MaterialLawBase):
    pass
