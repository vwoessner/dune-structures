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
        # TODO: Get rid if this non-sense!
        try:
            cell = u.ufl_element().cell()
        except:
            cell = u.ufl_operands[0].ufl_operands[0].ufl_element().cell()
        return [UFLPhysicalParameter(i, cell) for i in range(len(self.param_names))]

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
        dim = u.ufl_element().cell().topological_dimension()
        return Identity(dim) + grad(u)

    def infinitesimal_strain(self, u):
        return Variable(0.5 * (grad(u) + grad(u).T), label=Label(42))

    def right_cauchy_green(self, u):
        F = self.deformation_gradient(u)
        return dot(F.T, F)

    def left_cauchy_green(self, u):
        F = self.deformation_gradient(u)
        return dot(F, F.T)

    def cauchy_green_strain(self, u):
        C = self.right_cauchy_green(u)
        return Variable(0.5 * (C - Identity(3)), label=Label(43))

    def cauchy_green_invariants(self, u):
        B = self.left_cauchy_green(u)
        return [
            Variable(tr(B), label=Label(44)),
            Variable(0.5 * (tr(B) ** 2 - tr(dot(B, B))), label=Label(45)),
            Variable(det(B), label=Label(46)),
        ]

    def isochoric_cauchy_green_invariants(self, u):
        I1, I2, _ = self.cauchy_green_invariants(u)
        F = self.deformation_gradient(u)
        J = Variable(det(F), label=Label(47))
        return [
            Variable((J ** (-2.0 / 3.0)) * I1, label=Label(48)),
            Variable((J ** (-4.0 / 3.0)) * I2, label=Label(49)),
            J,
        ]


class InfinitesimalStrainBaseMaterialLaw(MaterialLawBase):
    def deformation_gradient(self, u):
        dim = u.ufl_element().cell().topological_dimension()
        return Identity(dim)

    def second_piola(self, u):
        return diff(self.strain_energy(u), self.infinitesimal_strain(u))


class CauchyGreenStrainBasedMaterialLaw(MaterialLawBase):
    def second_piola(self, u):
        return diff(self.strain_energy(u), self.cauchy_green_strain(u))


class CauchyGreenInvariantBasedMaterialLaw(MaterialLawBase):
    def second_piola(self, u):
        C = self.right_cauchy_green(u)
        I = self.cauchy_green_invariants(u)
        energy = self.strain_energy(u)

        gamma1 = diff(energy, I[0]) + I[0] * diff(energy, I[1])
        gamma2 = -diff(energy, I[1])
        gamma3 = I[2] * diff(energy, I[2])

        return 2.0 * (gamma1 * Identity(3) + gamma2 * C + gamma3 * inv(C))


class IsochoricCauchyGreenInvariantBasedMaterialLaw(MaterialLawBase):
    pass


class LinearMaterial(InfinitesimalStrainBaseMaterialLaw):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def id(self):
        return 0

    def strain_energy(self, u):
        epsilon = self.infinitesimal_strain(u)
        lmbda, mu = self.ufl_parameters(u)
        return 0.5 * lmbda * (tr(epsilon) ** 2) + mu * tr(epsilon * epsilon)


class StVenantKirchhoffMaterial(CauchyGreenStrainBasedMaterialLaw):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def id(self):
        return 1

    def strain_energy(self, u):
        E = self.cauchy_green_strain(u)
        lmbda, mu = self.ufl_parameters(u)
        return 0.5 * lmbda * (tr(E) ** 2) + mu * tr(dot(E, E))


class NeoHookeanMaterial(IsochoricCauchyGreenInvariantBasedMaterialLaw):
    @property
    def param_names(self):
        return ["first_lame", "second_lame"]

    @property
    def id(self):
        return 2

    def strain_energy(self, u):
        I1bar = self.isochoric_cauchy_green_invariants(u)[0]
        (param,) = self.ufl_parameters(u)
        return 0.5 * param * (I1bar - 3)

    # The following implementation is not yet based on isochoric invariants!
    # It should be removed once the base class IsochoricCauchyGreenInvariantBasedMaterialLaw
    # is operational...
    def cauchy_stress(self, u):
        lmbda, mu = self.ufl_parameters(u)
        F = self.deformation_gradient(u)
        B = self.left_cauchy_green(u)
        J = det(F)

        return lmbda * (J ** (-5.0 / 3.0)) * dev(B) + mu * (J - 1.0) * Identity(3)

    def first_piola(self, u):
        F = self.deformation_gradient(u)
        J = det(F)

        return J * self.cauchy_stress(u)


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
