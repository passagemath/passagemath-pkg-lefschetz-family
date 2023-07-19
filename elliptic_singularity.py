# -*- coding: utf-8 -*-

import sage.all

from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from sage.arith.misc import xgcd
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix

from sage.functions.other import floor
from sage.functions.other import ceil
from sage.arith.misc import gcd

from sage.symbolic.relation import solve
from sage.symbolic.expression import SR

from voronoi import FundamentalGroupVoronoi
from integrator import Integrator
from Util import Util
from Context import Context
from periods import LefschetzFamily

import logging
import time

logger = logging.getLogger(__name__)


class EllipticSingularities(object):
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

    fibre_confluence = {
        "I": [matrix(ZZ, [[1, 1], [0, 1]])],
        "II": [matrix(ZZ, [[2, 1], [-1, 0]]), matrix(ZZ, [[1, 1], [0, 1]])],
        "III": [matrix(ZZ, [[1, 0], [-1, 1]]), matrix(ZZ, [[2, 1], [-1, 0]]), matrix(ZZ, [[1, 1], [0, 1]])],
        "IV": [matrix(ZZ, [[2, 1], [-1, 0]]), matrix(ZZ, [[1, 1], [0, 1]]), matrix(ZZ, [[1, 0], [-1, 1]]), matrix(ZZ, [[2, 1], [-1, 0]])],
        "I*": [matrix(ZZ, [[-1,4],[-1,3]]), matrix(ZZ, [[1,0],[-1,1]])] + [matrix(ZZ, [[1,1],[0,1]])]*5,
        "II*": [matrix(ZZ, [[1, 0], [-1, 1]])]*8+[ matrix(ZZ, [[7, 1], [-36, -5]]), matrix(ZZ, [[4, 1], [-9, -2]])],
        "III*": [matrix(ZZ, [[1, 0], [-1, 1]])]*4 + [matrix(ZZ, [[4, 1], [-9, -2]])]*3 + [matrix(ZZ, [[3, 1], [-4, -1]])]*2,
        "IV*": [matrix(ZZ, [[1, 1], [0, 1]]), matrix(ZZ, [[1, 0], [-1, 1]]), matrix(ZZ, [[1, 1], [0, 1]]), matrix(ZZ, [[1, 0], [-1, 1]]), matrix(ZZ, [[2, 1], [-1, 0]]), matrix(ZZ, [[1, 1], [0, 1]]), matrix(ZZ, [[1, 1], [0, 1]]), matrix(ZZ, [[1, 0], [-1, 1]])],
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
        compl = (M-1).solve_right(vanishing*v)
        bc = matrix([vanishing, compl]).transpose()
        if bc.det()!=1:
            bc = matrix([-vanishing, compl]).transpose()
        assert bc.det()==1
        return bc, v
    
    @classmethod
    def normalize_Ivstar(cls, M):
        vanishing = (M+1).transpose().image().basis()[0]
        v = ZZ(gcd(vanishing))
        vanishing = vanishing/v
        vanishing = vanishing.change_ring(ZZ)
        compl = (M+1).solve_right(vanishing*v)
        bc = matrix([vanishing, compl]).transpose()
        if bc.det()!=1:
            bc = matrix([-vanishing, compl]).transpose()
        assert bc.det()==1
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
        possible_ds = range(ceil(solve(discr(d=x)>0,x)[0][0].right_hand_side()), floor(solve(discr(d=x)>0,x)[0][1].right_hand_side()))
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
        