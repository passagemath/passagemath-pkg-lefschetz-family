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
from sage.misc.flatten import flatten


from sage.symbolic.relation import solve
from sage.symbolic.ring import SR

from .monodromyRepresentation import MonodromyRepresentation
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

class MonodromyRepresentationFiberedProduct(MonodromyRepresentation):

    # def __init__(self, monodromy_matrices, intersection_product):
    #     super().__init__(monodromy_matrices, intersection_product)

    @property
    def types(self):
        if not hasattr(self, "_types"):
            types = []
            for M in self.monodromy_matrices:
                M1, M2 = self.disentangle(M)
                t, _, n = self.monodromy_class(M1)
                type1 = t + str(n) if t in ["I", "I*"] else t
                t, _, n = self.monodromy_class(M2)
                type2 = t + str(n) if t in ["I", "I*"] else t
                types += [type1 + " x " + type2]
            self._types = types
        return self._types
    
    def disentangle(self, M):
        M1 = M.submatrix(0,0,2,2)
        if M1 ==0:
            M1 = M.submatrix(0,2,2,2)
        M1 = (M1 / M1.det()**(1/2)).change_ring(ZZ)
        M2 = matrix(2, [(M.submatrix(i,j,2,2)*M1.inverse())[0,0] for i,j in [[0,0],[0,2],[2,0],[2,2]]])
        if M1==-identity_matrix(2) or M2==-identity_matrix(2):
            M1, M2 = -M1, -M2
        return M1, M2
    
    def desingularise_matrix(self, M):
        M1, M2 = self.disentangle(M)
        ty, base_change, nu = self.monodromy_class(M1)
        decomposition = self.fibre_confluence[ty]

        mats = []
        for M in decomposition[:-1]:
            mats += [tens11(base_change * M * base_change.inverse())]
            mats += [tens12(base_change * M * base_change.inverse())]
        mats += [tens11(base_change * decomposition[-1] * base_change.inverse()), tens12(base_change * decomposition[-1] * base_change.inverse())] * nu
        ty, base_change, nu = self.monodromy_class(M2)
        decomposition = self.fibre_confluence[ty]
        for M in decomposition[:-1]:
            mats += [tens21(base_change * M * base_change.inverse())]
            mats += [tens22(base_change * M * base_change.inverse())]
        mats += [tens21(base_change * decomposition[-1] * base_change.inverse()), tens22(base_change * decomposition[-1] * base_change.inverse())] * nu

        mats = [M.change_ring(ZZ) for M in mats]
        return mats
    
    def _compute_intersection_product_thimbles(self,i,j):
        vi = self.permuting_cycles_desingularisation[i]
        Mi = flatten(self.monodromy_matrices_desingularisation)[i]
        vj = self.permuting_cycles_desingularisation[j]
        Mj = flatten(self.monodromy_matrices_desingularisation)[j]
        di, dj = (Mi-1) * vi, (Mj-1) * vj
        res = di*self.fibre_intersection_product*dj
        resid = -vi*self.fibre_intersection_product*dj

        loopi = i//2
        loopj = j//2

        if loopi==loopj:
            return resid
        if i<j:
            return res
        else:
            return 0

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
    
    def lift(self, v):
            """Given a combination of thimbles of desingularisation, gives the corresponding homology class"""
            v = super(MonodromyRepresentationFiberedProduct, self).lift(v)
            vccts = [super(MonodromyRepresentationFiberedProduct, self).lift(v2) for v2 in flatten(self.vanishing_cycles_conifold_transition)]

            assert matrix(vccts) * self.intersection_product * v ==0, "cycle does not survive resolution"
        
            v = matrix(self.extensions_resolution + vccts).solve_left(v)
            
            return vector(list(v)[:-len(vccts)])
    
    @property
    def extensions_resolution(self):
        if not hasattr(self, "_extensions_resolution"):
            cosfs = matrix([super(MonodromyRepresentationFiberedProduct, self).lift(v) for v in flatten(self.vanishing_cycles_conifold_transition)])
            cosf_orth = (cosfs * self.intersection_product).change_ring(ZZ).right_kernel_matrix()
            toremove = cosfs.image().intersection(cosf_orth.image()).basis_matrix()
            compl = Util.find_complement(cosf_orth.solve_left(toremove).change_ring(ZZ))
            
            self._extensions_resolution = (compl * cosf_orth).rows()
        return self._extensions_resolution

    @property
    def homology(self): 
        if not hasattr(self, '_homology'):
            self._homology = identity_matrix(len(self.extensions_resolution)).rows()
        return self._homology

    @property
    def primary_lattice(self):
        if not hasattr(self, '_primary_lattice'):
            extensions = [super(MonodromyRepresentationFiberedProduct, self).lift(self.desingularise(v)) for v in self.extensions]
            components_of_singular_fibres = [super(MonodromyRepresentationFiberedProduct, self).lift(v) for v in flatten(self.components_of_singular_fibres, max_level=2)]
            primary_lattice = extensions + components_of_singular_fibres
            self._primary_lattice = matrix(primary_lattice).transpose()
        return self._primary_lattice


    @property
    def self_intersection_section(self):
        if not hasattr(self, '_self_intersection_section'):
            chi = ZZ((len(self.extensions_desingularisation)+ZZ(4))/ZZ(12))
            self._self_intersection_section = -chi
        return self._self_intersection_section
    
    @property
    def add(self):
        if not hasattr(self, '_add'):
            self._add = 0
        return self._add

    classes = {
        "I":matrix([[1,1],[0,1]]),
        "II":matrix([[1,1],[-1,0]]),
        "III":matrix([[0,1],[-1,0]]),
        "IV":matrix([[0,1],[-1,-1]]),
        "I*":matrix([[-1,-1],[0,-1]]),
        "II*":matrix([[0,-1],[1,1]]),
        "III*":matrix([[0,-1],[1,0]]),
        "IV*":matrix([[-1,-1],[1,0]]),
    }

    U = matrix([[1,1],[0,1]])
    V = matrix([[1,0],[-1,1]])

    fibre_confluence = {
        "I": [U],
        "II": [U, V],
        "III": [V,U,V],
        "IV": [U,V]*2,
        "I*": [U,V]*3+[U],
        "II*": [U,V]*5,
        "III*": [U,V]*3+[V,U,V],
        "IV*": [U,V]*4
    }

    @classmethod
    def monodromy_class(cls, M):
        M=M.change_ring(ZZ)
        if M.trace()==2:
            if M==identity_matrix(2):
                return "I", identity_matrix(2), 0
            else:
                return "I", cls.normalize_Iv(M)[0], cls.normalize_Iv(M)[1]
        if M.trace()==1:
            if M[0,1]>0:
                return "II", cls.find_conjugation(M, cls.classes["II"])[0], 1
            else:
                return "II*", cls.find_conjugation(M, cls.classes["II*"])[0], 1
        if M.trace()==0:
            try:
                bc = cls.find_conjugation(M, cls.classes["III"])[0]
                return "III", bc, 1
            except:
                bc = cls.find_conjugation(M, cls.classes["III*"])[0]
                return "III*", bc, 1
        if M.trace()==-1:
            try:
                bc = cls.find_conjugation(M, cls.classes["IV"])[0]
                return "IV", bc, 1
            except:
                bc = cls.find_conjugation(M, cls.classes["IV*"])[0]
                return "IV*", bc, 1
        if M.trace()==-2:
            if M==-identity_matrix(2):
                return "I*", identity_matrix(2), 0
            else:
                return "I*", cls.normalize_Ivstar(M)[0], cls.normalize_Ivstar(M)[1]
    

    @classmethod
    def normalize_Iv(cls, M):
        vanishing = (M-1).transpose().image().basis()[0]
        v = gcd(vanishing)
        v = ZZ(gcd(vanishing))
        vanishing = vanishing/v
        vanishing = vanishing.change_ring(ZZ)
        D, U, V = (M-1).smith_form()
        compl = V*vector([1,0])
        if (M-1)*compl != vanishing*v:
            compl = -compl
        assert (M-1)*compl == vanishing*v, "could not find permuting vector"
        bc = matrix([vanishing, compl]).transpose()
        if bc.det()!=1:
            bc = matrix([-vanishing, compl]).transpose()
        assert bc.det()==1, "could not find permuting vector"
        return bc, v
    
    @classmethod
    def normalize_Ivstar(cls, M):
        vanishing = (M+1).transpose().image().basis()[0]
        v = ZZ(gcd(vanishing))
        vanishing = vanishing/v
        vanishing = vanishing.change_ring(ZZ)
        D, U, V = (M+1).smith_form()
        compl = V*vector([1,0])
        if (M+1)*compl != vanishing*v:
            compl = -compl
        assert (M+1)*compl == vanishing*v, "could not find permuting vector"
        bc = matrix([vanishing, compl]).transpose()
        if bc.det()!=1:
            bc = matrix([-vanishing, compl]).transpose()
        assert bc.det()==1, "could not find permuting vector"
        return bc, v
    
    @classmethod
    def find_conjugation(cls, M1, M2):
        R = PolynomialRing(ZZ, ['a','b','c','d'])
        [a,b,c,d] = R.gens()
        bci = matrix([[a,b],[c,d]])
        eqs = (bci*M1 - M2*bci).coefficients() + [a*d-b*c-1]
        I = R.ideal(eqs)
        eq = I.elimination_ideal([a,b]).gens()[0]
        
        x = SR('x')
        discr = eq.coefficient(c)**2 - 4*eq.coefficient(c**2)*eq(c=0)
        possible_ds = range(ceil(solve(discr(d=x)>0,x)[0][0].right_hand_side()), floor(solve(discr(d=x)>0,x)[0][1].right_hand_side())+1)
        Qt = PolynomialRing(ZZ, 't')
        t=Qt.gens()[0]
        
        sols=[]
        for dv in possible_ds:
            possible_cs = eq(d=dv, c=t).roots(multiplicities=False)
            if len(possible_cs)>0:
                for cv in possible_cs:
                    I2 = R.ideal(eqs+[c-cv, d-dv])
                    av = -I2.elimination_ideal([b,c,d]).gens()[0](a=0)
                    bv = -I2.elimination_ideal([a,c,d]).gens()[0](b=0)
                    sols += [matrix([[av,bv],[cv,dv]])]    
        return [(s**-1).change_ring(ZZ) for s in sols]