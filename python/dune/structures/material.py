""" An abstraction for material laws for nonlinear mechanics.
This is inspired by a retired Simula project I came along:
https://bazaar.launchpad.net/~cbc-core/cbc.solve/main/files/head:/cbc/twist

New material models are defined in terms of a strain energy functional.
From that, the weak formulation is automatically derived.
"""

from dune.codegen.ufl.execution import *
from dune.structures.codegen import UFLPhysicalParameter
from ufl.classes import Variable, Label


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
        return Variable(0.5 * (grad(u) + grad(u).T), label=Label(42))

    def right_cauchy_green(self, u):
        F = self.deformation_gradient(u)
        return dot(F.T, F)

    def cauchy_green_strain(self, u):
        C = self.right_cauchy_green(u)
        return Variable(0.5 * (C - Identity(3)), label=Label(43))

    def cauchy_green_invariants(self, u):
        C = self.right_cauchy_green(u)
        return [Variable(tr(C), label=Label(44)),
                Variable(0.5*(tr(C)**2 - tr(dot(C, C))), label=Label(45)),
                Variable(det(C), label=Label(46))
                ]


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


class StVenantKirchhoffMaterial(MaterialLawBase):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def id(self):
        return 1

# Alternative - direct - formulation of 1st Piola-Kirchhoff stress tensor
#
#     def first_piola(self, u):
#         mu, lmbda = self.ufl_parameters(u)
#         F = self.deformation_gradient(u)
#         E = 0.5 * (dot(F, F.T) - Identity(3))
#         return 2.0 * mu * E + lmbda * tr(E) * Identity(3)
#

    def strain_energy(self, u):
        E = self.cauchy_green_strain(u)
        mu, lmbda = self.ufl_parameters(u)
        return 0.5 * lmbda * (tr(E)**2) + mu * tr(dot(E, E))

    def second_piola(self, u):
        return diff(self.strain_energy(u), self.cauchy_green_strain(u))


class CauchyGreenInvariantBasedMaterialLaw(MaterialLawBase):
    def second_piola(self, u):
        C = self.right_cauchy_green(u)
        I = self.cauchy_green_invariants(u)
        energy = self.strain_energy(u)

        gamma1 = diff(energy, I[0]) + I[0]*diff(energy, I[1])
        gamma2 = -diff(energy, I[1])
        gamma3 = I[2] * diff(energy, I[2])

        return 2.0 * (gamma1 * Identity(3) + gamma2 * C + gamma3 * inv(C))


class NeoHookeanMaterial(CauchyGreenInvariantBasedMaterialLaw):
    @property
    def param_names(self):
        return ["neoconst"]

    @property
    def id(self):
        return 2

    def strain_energy(self, u):
        I0 = self.cauchy_green_invariants(u)[0]
        param, = self.ufl_parameters(u)
        return param * (I0 - 3)


class MooneyRivlinMaterial(CauchyGreenInvariantBasedMaterialLaw):
    @property
    def param_names(self):
        return ["foo", "bar"]

    @property
    def id(self):
        return 3

    def strain_energy(self, u):
        I0, I1, _ = self.cauchy_green_invariants(u)
        c0, c1 = self.ufl_parameters(u)

        return c0 * (I0 - 3) + c1 * (I1 - 3)
