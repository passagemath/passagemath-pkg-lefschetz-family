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

from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import block_diagonal_matrix
from sage.matrix.special import zero_matrix

from sage.functions.other import floor
from sage.functions.other import ceil
from sage.arith.misc import gcd
from sage.misc.misc_c import prod
from sage.misc.flatten import flatten


from sage.symbolic.relation import solve
from sage.symbolic.ring import SR

from .monodromyRepresentation import MonodromyRepresentation
from .ellipticSingularity import EllipticSingularities
from .util import Util


import logging
import time

logger = logging.getLogger(__name__)

def tens1(M):
    return block_diagonal_matrix([M, M], subdivide=False)
def tens2(M):
    return block_matrix(list(map(lambda r:list(map(lambda c: c*identity_matrix(M.nrows()), r)), M.rows())), subdivide=False)

def tens11(M):
    return block_diagonal_matrix([M, identity_matrix(2)], subdivide=False)
def tens12(M):
    return block_diagonal_matrix([identity_matrix(2), M], subdivide=False)
def tens21(M):
    res = identity_matrix(4)
    for i in range(2):
        for j in range(2):
            res[2*i,2*j] = M[i,j]
    return res
def tens22(M):
    res = identity_matrix(4)
    for i in range(2):
        for j in range(2):
            res[2*i+1,2*j+1] = M[i,j]
    return res

class MonodromyRepresentationFibreProduct(MonodromyRepresentation):

    def __init__(self, monodromy_matrices, intersection_product, expected_types):
        super().__init__(monodromy_matrices, intersection_product)
        self.expected_types = expected_types

    @property
    def types(self):
        if not hasattr(self, "_types"):
            types = []
            for M, expected_type in zip(self.monodromy_matrices, self.expected_types):
                M1, M2 = self.disentangle(M)
                t, _, n = EllipticSingularities.monodromy_class(M1)
                if t + str(n) != expected_type:
                    M1, M2 = -M1, -M2
                    t, _, n = EllipticSingularities.monodromy_class(M1)
                assert t + str(n) == expected_type, "Unexpected fibre type"
                type1 = t + str(n) if t in ["I", "I*"] else t
                t, _, n = EllipticSingularities.monodromy_class(M2)
                type2 = t + str(n) if t in ["I", "I*"] else t
                types += [type1 + " x " + type2]
            self._types = types
        return self._types
    
    @staticmethod
    def disentangle(M):
        M1 = M.submatrix(0,0,2,2)
        if M1 ==0:
            M1 = M.submatrix(0,2,2,2)
        M1 = (M1 / M1.det()**(1/2)).change_ring(ZZ)
        M2 = matrix(2, [(M.submatrix(i,j,2,2)*M1.inverse())[0,0] for i,j in [[0,0],[0,2],[2,0],[2,2]]])
        return M1, M2
    
    @staticmethod
    def desingularise_matrix(M, expected_type):
        M1, M2 = MonodromyRepresentationFibreProduct.disentangle(M)
        ty, base_change, nu = EllipticSingularities.monodromy_class(M1)
        
        if ty + str(nu) != expected_type:
            M1, M2 = -M1, -M2
            ty, base_change, nu = EllipticSingularities.monodromy_class(M1)

        decomposition = EllipticSingularities.fibre_confluence[ty]

        mats = []
        for M in decomposition[:-1]:
            mats += [tens11(base_change * M * base_change.inverse())]
            mats += [tens12(base_change * M * base_change.inverse())]
        mats += [tens11(base_change * decomposition[-1] * base_change.inverse()), tens12(base_change * decomposition[-1] * base_change.inverse())] * nu
        ty, base_change, nu = EllipticSingularities.monodromy_class(M2)
        decomposition = EllipticSingularities.fibre_confluence[ty]
        for M in decomposition[:-1]:
            mats += [tens21(base_change * M * base_change.inverse())]
            mats += [tens22(base_change * M * base_change.inverse())]
        mats += [tens21(base_change * decomposition[-1] * base_change.inverse()), tens22(base_change * decomposition[-1] * base_change.inverse())] * nu

        mats = [M.change_ring(ZZ) for M in mats]
        return mats
    
    @property
    def monodromy_matrices_desingularisation(self):
        if not hasattr(self, '_monodromy_matrices_desingularisation'):
            monodromy_matrices_desingularisation = []
            for M, expected_type in zip(self.monodromy_matrices, self.expected_types):
                decomposition = self.desingularise_matrix(M, expected_type)
                assert prod(list(reversed(decomposition))) == M
                monodromy_matrices_desingularisation += [decomposition]
            self._monodromy_matrices_desingularisation = monodromy_matrices_desingularisation
        return self._monodromy_matrices_desingularisation


    def _compute_intersection_product_thimbles(self,i,j):
        vi = self.permuting_cycles_desingularisation[i]
        vj = self.permuting_cycles_desingularisation[j]

        Mi = flatten(self.monodromy_matrices_desingularisation)[i]
        Mj = flatten(self.monodromy_matrices_desingularisation)[j]
        
        di, dj = (Mi-1) * vi, (Mj-1) * vj

        res = di * self.fibre_intersection_product * dj
        resid = -vi * self.fibre_intersection_product * dj

        loopi = i//2
        loopj = j//2

        if loopi==loopj:
            return resid
        if loopi<loopj:
            return res
        else:
            return 0

    @property
    def intersection_product_resolution(self):
        if not hasattr(self, "_intersection_product_resolution"):
            cycles = matrix(self.extensions_resolution)
            self._intersection_product_resolution = cycles * self.intersection_product * cycles.transpose()
        return self._intersection_product_resolution

    @property
    def vanishing_cycles_conifold_transition(self):
        if not hasattr(self, "_vanishing_cycles_conifold_transition"):
            # first we figure out which pairs of singular components are there
            ords = []
            for Ms in self.monodromy_matrices_desingularisation:
                managed=False
                for i, M in enumerate(Ms):
                    if M.submatrix(2,0,2,2)!=0 or M.submatrix(0,2,2,2)!=0:
                        ords +=[i]
                        managed=True
                        break
                if not managed:
                    ords += [len(Ms)]
            # then we pair the vanishing cycles accordingly
            ranktot = 0
            rankmax = sum([len(Ms) for Ms in self.monodromy_matrices_desingularisation])
            sing_comps = []
            for n, M in zip(ords, self.vanishing_cycles_desingularisation):
                M = matrix(M)
                rank = M.dimensions()[0]
                vanishing = []
                Ni = M[:n].nrows()//2
                Nj = M[n:].nrows()//2
                for i in range(Ni):
                    for j in range(Nj):
                        vc = M[:n][2*i:2*i+2].stack(M[n:][2*j:2*j+2]).left_kernel_matrix().row(0)
                        vc = vector([0]*2*i + list(vc[:2]) + [0]*2*(Ni-i-1) + [0]*2*j +  list(vc[2:]) + [0]*2*(Nj-j-1))
                        vanishing += [vc]
                sing_comps += [[vector([0]*ranktot+list(v) + [0]*(rankmax-ranktot-rank)) for v in vanishing]]
                ranktot+=rank
            self._vanishing_cycles_conifold_transition = sing_comps
        return self._vanishing_cycles_conifold_transition
    
    def resolve(self, v):
            """Given a combination of thimbles of desingularisation, gives the corresponding homology class"""
            vccts = [self.lift(v2) for v2 in flatten(self.vanishing_cycles_conifold_transition)]
            assert matrix(vccts) * self.intersection_product * v ==0, "cycle does not survive resolution"
            v = matrix(self.extensions_resolution + vccts).solve_left(v)
            return vector(list(v)[:-len(vccts)])
    
    def proj(self, v):
            """Given a resolved homology class, gives a corresponding homology class (up to a conifold transition vanishing cycle)"""
            return v * matrix(self.extensions_resolution)
    
    @property
    def extensions_resolution(self):
        if not hasattr(self, "_extensions_resolution"):
            vanishing_cycles = matrix([self.lift(v) for v in flatten(self.vanishing_cycles_conifold_transition)])
            vanishing_cycles_orth = (vanishing_cycles * self.intersection_product).change_ring(ZZ).right_kernel_matrix()
            toremove = vanishing_cycles.image().intersection(vanishing_cycles_orth.image()).basis_matrix()
            compl = Util.find_complement(vanishing_cycles_orth.solve_left(toremove).change_ring(ZZ))
            
            self._extensions_resolution = (compl * vanishing_cycles_orth).rows()
        return self._extensions_resolution

    @property
    def homology(self): 
        if not hasattr(self, '_homology'):
            self._homology = identity_matrix(len(self.extensions_resolution)).rows()
        return self._homology

    @property
    def primary_lattice(self):
        if not hasattr(self, '_primary_lattice'):
            extensions = [super(MonodromyRepresentationFibreProduct, self).lift(self.desingularise(v)) for v in self.extensions]
            # extensions = matrix(extensions).change_ring(ZZ).image().basis()
            components_of_singular_fibres = [super(MonodromyRepresentationFibreProduct, self).lift(v) for v in flatten(self.components_of_singular_fibres, max_level=2)]
            primary_lattice = extensions + components_of_singular_fibres
            self._primary_lattice = matrix(primary_lattice).transpose()
        return self._primary_lattice
    
    @property
    def primary_periods(self):
        if not hasattr(self, '_primary_periods'):
            homology_mat = matrix(self.extensions).transpose()
            integrated_thimbles =  matrix(self.integrated_thimbles)
            self._primary_periods = integrated_thimbles * homology_mat
        return self._primary_periods


    
    @property
    def permuting_cycles(self):
        if not hasattr(self, '_permuting_cycles'):
            permuting_cycles = [identity_matrix(self.dim).rows() for i in range(len(self.monodromy_matrices))]
            self._permuting_cycles = permuting_cycles
        return self._permuting_cycles
    
    @property
    def infinity_loops(self):
        """The linear combinations of thimbles that correspond to extensions along the (trivial) loop around infinity."""
        if not hasattr(self, '_infinity_loops'):
            infinity_cycles = []
            for i in range(self.dim):
                v = vector([1 if k==i else 0 for k in range(self.dim)])
                coefs = []
                for j in range(len(self.monodromy_matrices)):
                    M = self.monodromy_matrices[j]
                    if len(self.permuting_cycles[j])==0:
                        continue
                    coefs += list(v)
                    v = self.monodromy_matrices[j]*v
                infinity_cycles+=[vector(coefs)]
            self._infinity_loops = matrix(infinity_cycles).change_ring(ZZ).rows()
        return self._infinity_loops
    
    @property
    def add(self):
        if not hasattr(self, '_add'):
            self._add = 0
        return self._add