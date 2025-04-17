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
from sage.functions.other import factorial
from sage.functions.other import floor
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import block_diagonal_matrix
from sage.arith.functions import lcm


from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.schemes.toric.weierstrass import WeierstrassForm
from sage.misc.flatten import flatten

from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice

from .voronoi import FundamentalGroupVoronoi
from .integrator import Integrator
from .util import Util
from .context import Context
from .monodromyRepresentation import MonodromyRepresentation
from sage.functions.other import binomial


import logging
import time

logger = logging.getLogger(__name__)

def tens1(M):
    return block_diagonal_matrix([M, M], subdivide=False)
def tens2(M):
    return block_matrix(list(map(lambda r:list(map(lambda c: c*identity_matrix(M.nrows()), r)), M.rows())), subdivide=False)
def tensor_prod(a,b):
    return vector([a[0]*b[0], a[1]*b[0], a[0]*b[1], a[1]*b[1]])

def derivatives_coords(S, w, k):
    derivatives = [S.P.parent()(0), w]
    for k in range(k):
        derivatives += [S._derivative(derivatives[-1], S.P)] 
    return S.family._coordinates(derivatives, S.basepoint)

def derivative_of_product(a, b, n):
    return sum([binomial(n, i) * tensor_prod(a.row(i+1), b.row(n-i+1)) for i in range(n+1)])

class FiberedProduct(object):
    def __init__(self, S1, S2, correction=ZZ(1), **kwds) -> None:
        """S1 and S2 are two elliptic surfaces with the same basepoint
        """
        
        self.ctx = Context(**kwds)
        
        self._S1 = S1
        self._S2 = S2

        self.correction = correction

        assert S1.basepoint == S2.basepoint, "the basepoint of the input elliptic surfaces should be the same"
        self._basepoint=S1.basepoint
    
    
    @property
    def S1(self):
        return self._S1
    @property
    def S2(self):
        return self._S2
    
    @property
    def basepoint(self):
        return self._basepoint

    @property
    def critical_values(self):
        if not hasattr(self, "_critical_values"):
            self._critical_values =[c for c in Util.remove_duplicates(self.S1.critical_values + self.S2.critical_values) if c != "infinity"]
        return self._critical_values

    @property
    def fundamental_group(self):
        if not hasattr(self,'_fundamental_group'):
            fundamental_group = FundamentalGroupVoronoi([e for e in self.critical_values if e!= 'infinity'], self.basepoint)
            fundamental_group.sort_loops()
            self._critical_values = fundamental_group.points[1:]
            self._fundamental_group = fundamental_group
        return self._fundamental_group
    
    @property
    def paths(self):
        if not hasattr(self, '_paths'):
            self._paths = [list(map(lambda i:self.fundamental_group.vertices[i], p)) for p in self.fundamental_group.pointed_loops]
        return self._paths

    @property
    def thimbles(self):
        return self.monodromy_representation.thimbles

    @property
    def permuting_cycles(self):
        return self.monodromy_representation.permuting_cycles
    
    @property
    def extensions(self):
        """Representants of the extensions of the elliptic surface."""
        return self.monodromy_representation.extensions

    @property
    def borders_of_thimbles(self):
        return self.monodromy_representation.borders_of_thimbles

    @property
    def infinity_loops(self):
        """The linear combinations of thimbles that correspond to extensions along the (trivial) loop around infinity."""
        return self.monodromy_representation.infinity_loops
        
    @property
    def vanishing_cycles(self):
        if not hasattr(self, "_vanishing_cycles"):
            self._vanishing_cycles = [[(M-1)*v for v in permuting_cycles_of_M] for permuting_cycles_of_M, M in zip(self.permuting_cycles, self.monodromy_matrices)]
        return self._vanishing_cycles

    @property
    def infinity_loops(self):
        """The linear combinations of thimbles that correspond to extensions along the (trivial) loop around infinity."""
        return self.monodromy_representation.infinity_loops

    @property
    def homology(self):
        return self.monodromy_representation.extensions
    
    def integrate(self, L):
        logger.info("Computing numerical transition matrices of operator of order %d and degree %d (%d edges total)."% (L.order(), L.degree(), len(self.fundamental_group.edges)))
        begin = time.time()

        integrator = Integrator(self.fundamental_group, L, self.ctx.nbits)
        transition_matrices = integrator.transition_matrices
        
        end = time.time()
        duration_str = time.strftime("%H:%M:%S",time.gmtime(end-begin))
        logger.info("Integration finished -- total time: %s."% (duration_str))

        return transition_matrices
    
    @property
    def monodromy_matrices(self):
        if not hasattr(self,'_monodromy_matrices'):
            logger.info("Computing monodromy matrices")
            L1 = self.S1.cyclic_picard_fuchs_equation
            L2 = self.S2.cyclic_picard_fuchs_equation
            Ltot = L1.symmetric_product(L2)

            Dt, = Ltot.parent().gens()
            
            transition_matrices = self.integrate(Ltot*Dt)
            
            w1, w2 = self.S1.P.parent()(1), self.S2.P.parent()(1)
            d1 = derivatives_coords(self.S1, w1, 3)
            d2 = derivatives_coords(self.S2, w2, 3)
            derivatives_tensor = matrix([vector([0]*4)] + [derivative_of_product(d1, d2, i) for i in range(4)])
            
            pM = tens1(self.S1.fibre.period_matrix)*tens2(self.S2.fibre.period_matrix)
            
            integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(5)])
            init_cond = (integration_correction * derivatives_tensor).submatrix(1,0)
            
            monodromy_matrices = [pM.inverse() * init_cond.inverse() * M.submatrix(1,1) * init_cond * pM for M in transition_matrices]
            monodromy_matrices = [M.change_ring(ZZ) for M in monodromy_matrices]


            Mtot = identity_matrix(4)
            for M in monodromy_matrices:
                Mtot = M * Mtot
            if Mtot != identity_matrix(4):
                self._critical_values = self.critical_values + ["infinity"]

                transition_matrix_infinity = 1
                for M in transition_matrices:
                    transition_matrix_infinity = M * transition_matrix_infinity
                transition_matrices += [transition_matrix_infinity.inverse()]
                
                monodromy_matrices += [(Mtot.inverse()).change_ring(ZZ)]
                
                pathtot=[]
                for path in self.paths:
                    pathtot = pathtot+path
                self._paths += [list(reversed(Util.simplify_path(pathtot)))]

            self._monodromy_matrices = monodromy_matrices
        return self._monodromy_matrices

    @property
    def monodromy_representation(self):
        if not hasattr(self,'_monodromy_representation'):
            fibre_intersection_product = tens1(self.S1.fibre.intersection_product)*tens2(self.S2.fibre.intersection_product)
            self._monodromy_representation = MonodromyRepresentation(self.monodromy_matrices, fibre_intersection_product)
        return self._monodromy_representation

    @property
    def cohomology(self):
        """a cohomology class is represented by `(w1, w2, denom)` where `w1` is a form of `self.S1`, `w2` is a form of `self.S2` and `denom` is a rational function of the parameter."""
        if not hasattr(self,"_cohomology"):
            S = self.S1.P.parent()
            R = S.base_ring()
            
            Swithu = PolynomialRing(R, ['t', 'u'])
            t,u = Swithu.gens() 
            phi = S.hom([t*u])
            phii = Swithu.hom([S.gens()[0], 1])
            f1 = Swithu(1)/phi(self.S1.P) # f1 is the rational function of which the residue is the holomorphic form of S1
            
            derivatives = [f1]
            for i in range(3):
                derivatives += [derivatives[-1].derivative(u)]
            nums = [phii(phi(self.S1.P**(i+1))*w)/factorial(i) for i, w in enumerate(derivatives)]

            w2 = vector([1,0])
            w1s = [S(Swithu(n)) for n in nums]
            w1s, denom = self.S1.family.coordinates(w1s)

            self._cohomology = [(w1, w2, denom) for w1 in w1s]
        return self._cohomology

    @property
    def picard_fuchs_equations(self):
        if not hasattr(self,"_picard_fuchs_operators"):
            picard_fuchs_equations = []
            for w1, w2, denom in self.cohomology:
                L1 = self.S1.family.picard_fuchs_equation(w1)
                L2 = self.S2.family.picard_fuchs_equation(w2)
                L = L1.symmetric_product(L2)*denom
                L = L*(1/self.correction)
                mul = lcm([c.denominator() for c in L.coefficients()])
                L = mul*L
                L = L.change_ring(L.base_ring().ring_of_integers())
                picard_fuchs_equations += [L]
            self._picard_fuchs_equations = picard_fuchs_equations
        return self._picard_fuchs_equations
    
    @property
    def transition_matrices(self):
        if not hasattr(self, '_transition_matrices'):
            transition_matrices = []
            for L in self.picard_fuchs_equations:
                L = L* L.parent().gens()[0]
                transition_matrices += [self.integrate(L)]
                if "infinity" in self.critical_values:
                    transition_matrix_infinity = 1
                    for M in transition_matrices[-1]:
                        transition_matrix_infinity = M*transition_matrix_infinity
                    transition_matrices[-1] += [transition_matrix_infinity**(-1)]
            self._transition_matrices = transition_matrices
        return self._transition_matrices
    
    @property
    def integrated_thimbles(self):
        if not hasattr(self, '_integrated_thimbles'):
            integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(5)])
            pM = tens1(self.S1.fibre.period_matrix)*tens2(self.S2.fibre.period_matrix)

            _integrated_thimbles_all = []

            for transition_matrices, w in zip(self.transition_matrices, self.cohomology):
                w1, w2, denom = w

                numerator1 = (w1[0]*self.S1.P + w1[1]*self.S1.family.pol.parent()(self.S1.family.coho1.basis()[1]))*self.correction/denom
                numerator2 = (w2[0]*self.S2.P + w2[1]*self.S2.family.pol.parent()(self.S2.family.coho1.basis()[1]))
                d1 = derivatives_coords(self.S1, numerator1, 3)
                d2 = derivatives_coords(self.S2, numerator2, 3)
                derivatives_tensor = matrix([vector([0]*4)] + [derivative_of_product(d1, d2, i) for i in range(4)])
                initial_conditions = integration_correction * derivatives_tensor * pM

                _integrated_thimbles = []
                for i, ps in enumerate(self.permuting_cycles):
                    for p in ps:
                        _integrated_thimbles += [(transition_matrices[i] * initial_conditions * p)[0]]
                _integrated_thimbles_all += [_integrated_thimbles]
            self._integrated_thimbles = _integrated_thimbles_all
        return self._integrated_thimbles

    @property
    def period_matrix(self):
        if not hasattr(self, "_period_matrix"):
            homology_mat = matrix(self.extensions).transpose()
            integrated_thimbles =  matrix(self.integrated_thimbles)
            self._period_matrix = integrated_thimbles * homology_mat
        return self._period_matrix
    
    @property
    def intersection_product(self):
        return self.monodromy_representation.intersection_product_extensions
    
    # @property
    # def intersection_product(self):
    #     if not hasattr(self,'_intersection_product'):
    #         # in the next line, 12 is specific to this case
    #         inter_prod_thimbles = matrix([[self._compute_IP_thimbles(i,j) for j in range(12)] for i in range(12)])
    #         intersection_product = self.homology*inter_prod_thimbles*self.homology.transpose()
    #         intersection_product = intersection_product.change_ring(ZZ)

    #         self._intersection_product = intersection_product
    #     return self._intersection_product

    # def _compute_IP_thimbles(self, i, j):
    #     loops_ids = flatten([[i]*len(vanishing_cycles_of_M) for i, vanishing_cycles_of_M in enumerate(self.vanishing_cycles)])
    #     tensIP = tens1(self.S1.fibre.intersection_product)*tens2(self.S2.fibre.intersection_product)

    #     loopi, loopj = loops_ids[i], loops_ids[j]
    #     vi = flatten(self.permuting_cycles)[i]
    #     Mi = self.monodromy_matrices[loopi]
    #     vj = flatten(self.permuting_cycles)[j]
    #     Mj = self.monodromy_matrices[loopj]

    #     di, dj = (Mi-1)*vi, (Mj-1)*vj

    #     res = di*tensIP*dj
    #     resid = -vi*tensIP*dj

    #     if loopi == loopj:
    #         return resid
    #     if loopi < loopj:
    #         return res
    #     else:
    #         return 0