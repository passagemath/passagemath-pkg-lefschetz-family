# -*- coding: utf-8 -*-

# lefschetz-family
# Copyright (C) 2021  Eric Pichon-Pharabod

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sage.all

from .numperiods.family import Family
from .numperiods.integerRelations import IntegerRelations
from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.qqbar import QQbar
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import NumberField
from sage.rings.imaginary_unit import I
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.functions.other import factorial
from sage.functions.other import floor
from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import block_diagonal_matrix
from sage.matrix.special import zero_matrix
from sage.arith.functions import lcm
from sage.arith.misc import gcd
from sage.arith.misc import factor
from sage.symbolic.constants import pi
from sage.functions.transcendental import zeta
from sage.misc.misc_c import prod
from sage.groups.matrix_gps.symplectic import Sp

from sage.misc.flatten import flatten
from sage.misc.functional import log

from .voronoi import FundamentalGroupVoronoi
from .integrator import Integrator
from .util import Util
from .context import Context
from .hypersurface import Hypersurface
from .monodromyRepresentationEllipticSurface import MonodromyRepresentationEllipticSurface

import logging
import time

logger = logging.getLogger(__name__)


class CalabiYauOperator(object):
    def __init__(self, L, basepoint=None, **kwds):
        """L, a differential equation 
        """
        
        self.ctx = Context(**kwds)

        self._L = L
        
        if basepoint!= None: # it is useful to be able to specify the basepoint to avoid being stuck in arithmetic computations if critical values have very large modulus
            assert basepoint not in self.singular_values, "basepoint is not regular"
            self._basepoint = basepoint

        if not self.ctx.debug:
            fg = self.fundamental_group # this allows reordering the critical points straight away and prevents shenanigans. There should be a better way to do this
    
    @property
    def L(self):
        """returns the Calabi-Yau operator"""
        return self._L
    
    @property
    def order(self):
        return self.L.order()
    
    @property
    def fundamental_group(self):
        if not hasattr(self, "_fundamental_group"):
            fundamental_group = FundamentalGroupVoronoi(self.singular_values, self.basepoint)
            self._singular_values = [self._singular_values[i] for i in fundamental_group.sort_loops()]
            self._fundamental_group = fundamental_group
        return self._fundamental_group

    @property
    def loops(self):
        return [[self.fundamental_group.vertices[i] for i in loop] for loop in self.fundamental_group.pointed_loops]

    @property
    def singular_values(self):
        if not hasattr(self, "_singular_values"):
            self._singular_values = self.L.leading_coefficient().roots(QQbar, multiplicities=False)
        return self._singular_values

    @property
    def basepoint(self):
        if not hasattr(self, "_basepoint"):
            self._basepoint = QQ(10)**(min([floor(log(abs(r),10)) for r in self.singular_values if r!=0])-1)
        return self._basepoint

    @property
    def period_matrix(self):
        if not hasattr(self, "_period_matrix"):
            for starting_vector in ((self.transition_matrices[self.maximal_unipotent_monodromy_index]-1)**3).columns():
                if starting_vector!=0:
                    break
            Pi = self._discover_rational_basis(starting_vector)
            if not self.is_semi_simple:
                raise Exception("holomorphic period generates stable subsystem - operator is factorizable")
            try:
                monodromy_on_homology = self.monodromy_from_periods(Pi)
            except:
                raise NotImplementedError("Non-integral monodromy -- Look for algebraic coefficients")
            CB = Util.saturate(monodromy_on_homology) if self.has_integral_monodromy else identity_matrix(4)
            self._period_matrix =  Pi * CB
        return self._period_matrix
    
    @property
    def number_field(self):
        if not hasattr(self, "_number_field"):
            for starting_vector in ((self.transition_matrices[self.maximal_unipotent_monodromy_index]-1)**3).columns():
                if starting_vector!=0:
                    break
            Pi = self._discover_rational_basis(starting_vector)
            t = self.L.base_ring().gen(0)
            coefs = flatten([(Pi.inverse() * M * Pi).coefficients() for M in self.transition_matrices])
            unknowns = [1]
            for c in coefs:
                if IntegerRelations(matrix(unknowns+[c]).transpose()).basis.nrows()==0:
                    unknowns += [c]

            self._has_integral_monodromy = len(unknowns) == 1

            if not self.has_integral_monodromy:
                P = lcm([sum([c*t**i for i, c in enumerate(Util.check_if_algebraic(u))]) for u in unknowns])
                SF = P.splitting_field("nu").optimized_representation(name="nu")[0]
                P = SF.defining_polynomial()

                r = P.roots(QQbar, multiplicities=False)[1]
                SF = NumberField(P,"nu", embedding=r)
                self._number_field = (SF, r)
            else: 
                self._number_field = (QQ, 1)
        return self._number_field
    
    def _discover_rational_basis(self, starting_vector):
        vs = [starting_vector]
        while len(vs)!=4:
            found = False
            for M in self.transition_matrices:
                for v in vs:
                    Pi = matrix(vs + [M*v]).transpose()
                    if any([m!=0 for m in Pi.minors(Pi.ncols())]):
                        vs += [M*v]
                        found = True
            if not found:
                self._is_semi_simple = False
                break
        if found:
            self._is_semi_simple = True
        return matrix(vs).transpose()

    @property
    def intersection_product(self):
        return self.compute_intersection_product(self.monodromy_matrices)

    def compute_intersection_product(self, Ms):
        SF, _ = self.number_field
        R = PolynomialRing(SF, self.order**2, 'a')
        IP = matrix(self.order,self.order,R.gens())
        I = R.ideal(flatten((IP+IP.transpose()).coefficients()+[(M.transpose()*IP*M-IP).coefficients() for M in Ms]))
        gens = I.groebner_basis()
        for a in R.gens():
            if a not in gens:
                break
        I2 = R.ideal(flatten([a-1] + (IP+IP.transpose()).coefficients()+[(M.transpose()*IP*M-IP).coefficients() for M in Ms]))
        IP = IP.apply_map(lambda c: I2.reduce(c)).change_ring(self.number_field[0])
        IP = IP * lcm([c.denominator() for c in IP.coefficients()])
        if self.has_integral_monodromy:
            IP = IP.change_ring(ZZ)
        return IP


    @property
    def maximal_unipotent_monodromy_index(self):
        """returns the Calabi-Yau operator"""
        return self._fundamental_group.points.index(0)-1 # TODO : what if we want a different MUM_point ?

    @property
    def transition_matrices(self):
        if not hasattr(self, "_transition_matrices"):
            integrator = Integrator(self.fundamental_group, self.L, self.ctx.nbits)
            transition_matrices = integrator.transition_matrices
            if self.L.annihilator_of_composition(1/self.L.base_ring().gen()).leading_coefficient()(0)==0:
                transition_matrices += [prod(list(reversed(transition_matrices))).inverse()]
                self._singular_values = self.singular_values + ["infinity"]
            self._transition_matrices = transition_matrices
        return self._transition_matrices
    
    @property
    def types(self):
        if not hasattr(self, "_types"):
            types = [self._monodromy_type(M) for M in self.transition_matrices]
            self._types = types
        return self._types
    
    @classmethod
    def _monodromy_type(cls, M):
        for i in range(1,13):
            for j in range(1,5):
                try:
                    if ((M**i-1)**j).change_ring(ZZ) == 0:
                        return (i,j)
                except KeyboardInterrupt:
                    raise KeyboardInterrupt
                except:
                    continue
        raise Exception("Monodromy type not identified")

    @property
    def maximal_unipotent_monodromy_fibres(self):
        if not hasattr(self, "_maximal_unipotent_monodromy_fibres"):
            self._maximal_unipotent_monodromy_fibres = [i for i, t in enumerate(self.types) if t == (1,4)]
        return self._maximal_unipotent_monodromy_fibres
    
    @property
    def conifold_fibres(self):
        if not hasattr(self, "_conifold_fibres"):
            conifold_fibres = [i for i, t in enumerate(self.types) if t == (1,2)]
            conifold_fibres.sort(key = lambda i: abs(self.singular_values[i])if self.singular_values[i]!="infinity" else max([abs(c) for c in self.singular_values if c!="infinity"])+1)
            self._conifold_fibres = conifold_fibres
        return self._conifold_fibres

    def _reduce_discriminant(self, vs, monodromy_matrices):
        CB = matrix(vs).transpose()
        discr = self.compute_intersection_product([CB.inverse() * M * CB for M in monodromy_matrices]).det()
        managed=True
        while managed:
            managed=False
            for p in discr.divisors():
                for v in vs:
                    v2s = [v2 if v2!=v else v/p for v2 in vs]
                    discr2 = self._extend(v2s, monodromy_matrices)
                    if discr2<discr:
                        vs = v2s
                        discr = discr2
                        managed = True
                        break
                if managed:
                    break
        return vs

    def _extend(self, vs, Ms):
        CB = matrix(vs).transpose()
        newMs = [(CB.inverse()*M*CB) for M in Ms]
        if self.has_integral_monodromy:
            newCB = CB*Util.saturate(newMs)
            newMs = [(newCB.inverse()*M*newCB).change_ring(ZZ) for M in Ms]
        IP = self.compute_intersection_product(newMs).change_ring(ZZ)
        return IP.det()

    @property
    def basis_change_to_mum_point(self):
        if not hasattr(self, "_basis_change_to_mum_point"): # performance can be gained here (redundant computation with self.transition_matrices)
            path = [self.fundamental_group.vertices[i] for i in self.fundamental_group.paths[self.maximal_unipotent_monodromy_index]] + [0]
            self._basis_change_to_mum_point = self.L.numerical_transition_matrix(path, ZZ(2)**-self.ctx.nbits, assume_analytic=True)
        return self._basis_change_to_mum_point

    @property
    def gamma_class(self):
        if not hasattr(self, "_gamma_class"):
            gamma_class = self.compute_decomposition_in_mum_frobenius(self.period_matrix * self.standard_basis)
            self._gamma_class = gamma_class
        return self._gamma_class
    
    @property
    def monodromy_matrices(self):
        if not hasattr(self, "_monodromy_matrices"):
            standard_periods = self.period_matrix * self.standard_basis
            self._monodromy_matrices = self.monodromy_from_periods(standard_periods)
        return self._monodromy_matrices

    @property
    def laddered_basis(self):
        if not hasattr(self, "_laddered_basis"):
            monodromy_matrices = self.monodromy_from_periods(self.period_matrix)
            Mmum = monodromy_matrices[self.maximal_unipotent_monodromy_index]
            holo = ((Mmum-1)).right_kernel().gen(0)

            found = False
            for i in self.conifold_fibres:
                M = monodromy_matrices[i]
                if M*holo not in ((Mmum-1)**3).right_kernel():
                    log4 = M*holo
                    found = True
                    break

            if not found:
                for M in monodromy_matrices:
                    if M*holo not in ((Mmum-1)**3).right_kernel():
                        log4 = M*holo
                        found = True
                        break
            if found:
                d = (2*3*5*7*11*13*17*19*21*23*29)**10
                holo, log2, log3 = (Mmum-1)**3*log4/d, (Mmum-1)**2*log4/d, (Mmum-1)*log4

            if not found:
                vs = [holo]
                log2ker =((Mmum-1)**2).right_kernel_matrix()
                log2 = vector(Util.find_complement(log2ker.solve_left(matrix(vs)), primitive=False) * log2ker)
                vs += [log2]
                
                log3ker =((Mmum-1)**3).right_kernel_matrix()
                log3 = vector(Util.find_complement(log3ker.solve_left(matrix(vs)), primitive=False) * log3ker)
                vs += [log3]
                
                log4 = vector(Util.find_complement(matrix(vs), primitive=False))
                vs += [log4]

            vs = [log4, log3, log2, holo]
            vs = matrix(vs).transpose()
            
            CB = Util.saturate([vs.inverse() * M * vs for M in monodromy_matrices]) if self.has_integral_monodromy else identity_matrix(4)
            perm = matrix(ZZ, 4, 4, {(0,0):1, (1,1):1, (2,3):1, (3,2):1})
            self._laddered_basis = (vs * CB * perm)
        return self._laddered_basis
    

    def normalize_basis(self, basis):
        monodromy_matrices = self.monodromy_from_periods(self.period_matrix*basis)
        if self.has_integral_monodromy:
            IP = self.compute_intersection_product(monodromy_matrices)
            vs = identity_matrix(4).rows()
            vs[3] -= vs[2] * ((vs[0] * IP * vs[3])//(vs[0] * IP * vs[2]))
            vs[1] -= vs[2] * ((vs[0] * IP * vs[1])//(vs[0] * IP * vs[2]))
            vs = matrix(vs).transpose()

            res = self.compute_decomposition_in_mum_frobenius(self.period_matrix * basis * vs)

            vs = vs.columns()
            vs[3] -= ((QQ(res[3,0]).numerator() * QQ(res[2,0]).denominator()) // (QQ(res[2,0]).numerator() * QQ(res[3,0]).denominator()))*vs[2]
            vs[0] -= vs[1] * ((vs[0] * IP * vs[3])//(vs[1] * IP * vs[3]))
            vs[3] -= vs[2] * ((vs[0] * IP * vs[3])//(vs[0] * IP * vs[2]))

            if res[3,1]<0:
                vs[3] = -vs[3]
            if vs[3] * IP * vs[1] >0:
                vs[1] = -vs[1]
            if res[2,0]<0:
                vs[2] = -vs[2]
            if vs[2] * IP * vs[0] >0:
                vs[0] = -vs[0]
            
            vs = matrix(vs).transpose()

            res = self.compute_decomposition_in_mum_frobenius(self.period_matrix * basis * vs)

            vs = vs.columns()
            vs[0] -= floor(res[0,0].monomial_coefficient(res[0,0].parent()(1))) * vs[2]
            vs[1] -= (QQ(res[1,1]).numerator() * QQ(res[3,1]).denominator()) // (QQ(res[3,1]).numerator() * QQ(res[1,1]).denominator()) * vs[3]
            vs = matrix(vs).transpose()

            return basis * vs
        return basis
    
    def monodromy_from_periods(self, periods):
        SF, r = self.number_field
        nu, d = SF.gen(0), SF.degree()
        variables = [r**i for i in range(d)]
        constants = [nu**i for i in range(d)]
        return [(periods.inverse() * M * periods).apply_map(lambda c: Util.get_coefficient(c, variables, constants)) for M in self.transition_matrices]

    @property
    def standard_basis(self):
        perm = matrix(ZZ, 4, 4, {(0,0):1, (1,1):1, (2,3):1, (3,2):1})
        basis = self.normalize_basis(self.laddered_basis) * perm

        monodromy_matrices = self.monodromy_from_periods(self.period_matrix * basis)
        CB = identity_matrix(4)
        if self.has_integral_monodromy: #the following is outdated 
            # vs = self._reduce_discriminant(identity_matrix(4).rows(), monodromy_matrices)
            vs = identity_matrix(4)
            CB = matrix(vs).transpose()
            CB = CB * Util.saturate([CB.inverse() * M * CB for M in monodromy_matrices])

        basis = self.normalize_basis(basis * CB * perm)

        basis = basis * lcm([c.denominator() for c in basis.coefficients()])

        return basis


    def compute_decomposition_in_mum_frobenius(self, periods):
        twopiI = self.ctx.CBF(ZZ(2)*pi*I)
        finalCB = (diagonal_matrix([1, twopiI**-1, twopiI**-2, twopiI**-3]) * matrix(ZZ, 4, 4, {(0,3):1, (1,2):1, (2,1):1, (3,0):1}) * diagonal_matrix([6,2,1,1])).transpose()
        period_matrix = self.basis_change_to_mum_point * periods

        Mmum = (self.monodromy_from_periods(periods)[self.maximal_unipotent_monodromy_index]-1)
        if self.has_integral_monodromy:
            Mmum = Mmum.change_ring(ZZ)
        c = period_matrix * Mmum.right_kernel().gen(0)
        period_matrix = period_matrix / c[-1]
        M = (finalCB.inverse()*period_matrix).transpose()

        SF, r = self.number_field
        nu = SF.gen(0)
        R = PolynomialRing(SF, 'lambda')
        constants = [1, self.ctx.CBF(zeta(3))/twopiI**3]
        variables = [1, R.gen(0)]

        constantsall, variablesall = [], []
        for i in range(SF.degree()):
            constantsall += [self.ctx.CBF(r)**i * v for v in constants]
            variablesall += [nu**i * v for v in variables]

        formal_matrix = M.apply_map(lambda c: Util.get_coefficient(c, constantsall, variablesall))

        return formal_matrix
    

    @property
    def monodromy_in_scaled_frobenius_basis(self):
        twopiI = self.ctx.CBF(ZZ(2)*pi*I)
        scaled_frobenius_basis = (diagonal_matrix([1, twopiI**-1, twopiI**-2, twopiI**-3]) * matrix(ZZ, 4, 4, {(0,3):1, (1,2):1, (2,1):1, (3,0):1}) * diagonal_matrix([6,2,1,1])).transpose()
        period_matrix = self.basis_change_to_mum_point.inverse() * scaled_frobenius_basis

        monodromy_matrices = [period_matrix.inverse() * M * period_matrix for M in self.transition_matrices]
        
        SF, r = self.number_field
        nu = SF.gen(0)
        R = PolynomialRing(SF, 'lambda')
        l = R.gen(0)
        variables = [1, l, l**2]
        hom = R.hom([self.ctx.CBF(zeta(3))/twopiI**3])
        constants = [hom(v) for v in variables]

        constantsall, variablesall = [], []
        for i in range(SF.degree()):
            constantsall += [self.ctx.CBF(r)**i * v for v in constants]
            variablesall += [nu**i * v for v in variables]

        monodromy_matrices = [M.apply_map(lambda c: Util.get_coefficient(c, constantsall, variablesall)) for M in monodromy_matrices]

        return monodromy_matrices
    

    @property
    def is_semi_simple(self):
        """returns whether the monodromy group is semi-simple"""
        if not hasattr(self, "_is_semi_simple"):
            for starting_vector in ((self.transition_matrices[self.maximal_unipotent_monodromy_index]-1)**3).columns():
                if starting_vector!=0:
                    break
            self._discover_rational_basis(starting_vector)
        return self._is_semi_simple
    
    @property
    def has_conifold(self):
        if not hasattr(self, "_has_conifold"):
            self._has_conifold = len(self.conifold_fibres)>0
        return self._has_conifold
    
    @property
    def is_symplectic(self):
        if not hasattr(self, "_is_symplectic"):
            self._is_symplectic = self.intersection_product == matrix(ZZ, 4, 4, {(0,2):1, (1,3):1, (2,0):-1, (3,1):-1})
        return self._is_symplectic
    
    @property
    def has_two_maximal_unipotent_monodromy_points(self):
        if not hasattr(self, "_has_two_maximal_unipotent_monodromy_points"):
            self._has_two_maximal_unipotent_monodromy_points = len(self.maximal_unipotent_monodromy_fibres)>1
        return self._has_two_maximal_unipotent_monodromy_points
    
    @property
    def has_integral_monodromy(self):
        if not hasattr(self, "_has_integral_monodromy"):
            self.number_field
        return self._has_integral_monodromy

    @property
    def canonical_topological_invariants(self):
        """ returns (H3, c2H, chi, p) where H3 is the triple intersection number, c2H is the second Chern class, chi is the top Chern class (the Euler characteristic), """
        raise NotImplemented

    @property
    def c2H(self):
        """returns the second Chern class of the operator"""
        raise NotImplemented

    @property
    def triple_intersection(self):
        """returns the triple intersection number of the operator"""
        raise NotImplemented

    @property
    def paths(self):
        if not hasattr(self, "_paths"):
            paths = [[self.fundamental_group.vertices[i] for i in path] for path in self.fundamental_group.pointed_loops]
            self.transition_matrices
            if "infinity" in self.singular_values:
                pathtot = []
                for path in paths:
                    pathtot+=path
                paths += [Util.simplify_path(pathtot)]
            self._paths = paths
        return self._paths
    
    def monodromy_index_mod_n(self, n):
        if not hasattr(self, "_reduction_indices"):
            self._reduction_indices = {}
        res = 1
        for p, power in factor(n):
            q = p**power
            if q not in self._reduction_indices.keys():
                self._reduction_indices[q] = self._monodromy_index_mod_q(q)
            res *= self._reduction_indices[q]
        return res

    def _monodromy_index_mod_q(self, q):
        R = IntegerModRing(q)
        monodromy_matrices_mod_q = [M.change_ring(R) for M in self.monodromy_matrices]
        Sp4 = Sp(4, R)
        if Sp4.invariant_form() == matrix(R, 4, 4, {(0,3):1, (1,2):1, (2,1):-1, (3,0):-1}):
            CB = matrix(R, 4, 4, {(0,0):1, (1,1):1, (2,3):1, (3,2):1})
            monodromy_matrices_mod_q = [CB.inverse()*M*CB for M in monodromy_matrices_mod_q]
        monodromy_group = Sp4.subgroup(monodromy_matrices_mod_q)
        return Sp4.cardinality()/monodromy_group.cardinality()