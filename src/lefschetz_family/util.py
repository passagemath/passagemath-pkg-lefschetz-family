# -*- coding: utf-8 -*-

import sage.all

from ore_algebra import *
from sage.rings.complex_mpfr import ComplexField
from sage.functions.other import floor
from sage.arith.misc import gcd
from sage.arith.misc import xgcd
from sage.combinat.integer_vector import IntegerVectors
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import identity_matrix

from sage.misc.prandom import randint


import logging

logger = logging.getLogger(__name__)


class Util(object):

    @classmethod 
    def simple_rational(cls, p, r):
        """Gives a rational q such that |p-q|<=r"""
        x=p
        l=[floor(x)]
        while abs(p-cls.evaluate_continued_fraction(l))>r:
            x=1/(x-l[-1])
            l+=[floor(x)]
        return cls.evaluate_continued_fraction(l)
    
    @classmethod
    def evaluate_continued_fraction(cls, l):
        """ Given a list l, evaluates the continued fraction l[0] + 1/(l[1] + 1/(l[2] + ...))"""
        p=l[-1]
        l=l[:-1]
        while len(l)>0:
            p= l[-1] +1/p
            l=l[:-1]
        return p

    @classmethod
    def invert_permutation(cls, l):
        """Given a list representing a permutation of [0, ..., len(l)-1], l[i] = j, returns the inverse permutation l2[j] = i"""
        return [l.index(x) for x in range(len(l))]

    @classmethod
    def simplify_path(cls, p):
        """Given a list of numbers p, returns a ``simplification'' of p, removing all the backtracking
        
        For example, `simplify_path([1,2,3,2,1,4,5,1]) == [1,4,5,1]`
        """
        i=1
        res = list(p)
        while i<len(res)-1:
            p = res[i-1]
            a = res[i]
            n = res[i+1] 
            if p==a:
                res = res[:i]+res[i+1:]
                if i!=1:
                    i=i-1
            elif p==n:
                res = res[:i]+res[i+2:]
                if i!=1:
                    i=i-1
            else:
                i=i+1
        return res

    @classmethod
    def monomials(cls, ring, degree):
        """Given a polynomial ring and an integer d, returns all the monomials of the ring with degree d."""
        return [ring.monomial(*m) for m in list(IntegerVectors(degree, ring.ngens()))]

    @classmethod
    def xgcd_list(cls, l):
        """Given a list of integers l, return a double consisting of the gcd of these numbers and the coefficients of a Bezout relation."""
        if len(l)==0:
            return 0, []
        if len(l)==1:
            return l[0], [1]
        d = gcd(l)
        result = [1]
        a = l[0]
        for i in range(len(l)-1):
            b = l[i+1]
            d2, u, v = xgcd(a,b)
            result = [k*u for k in result] + [v]
            a=d2
        assert d2==d, "not getting the correct gcd"
        return d2, result


    @classmethod
    def path(cls, path, x):
        """Given a list of complex numbers path, and a number 0<x<1, return the path(x) where path is seen as a path [0,1]\\to C."""
        CC=ComplexField(500)
        dtot = sum([CC(abs(p1-p2)) for (p1,p2) in zip(path[:-1], path[1:])])
        dmin, dmax = 0, CC(abs(path[0]-path[1]))
        for i in range(len(path)-1):
            if x*dtot<=dmax and x*dtot>=dmin:
                break;
            else:
                dmin, dmax=dmax, dmax+CC(abs(path[i+1]-path[i+2]))
        t = Util.simple_rational((x*dtot -dmin)/(dmax-dmin), 10e-10)
        return (1-t)*path[i] + t*path[i+1]



    @classmethod
    def select_closest(cls, l, e):
        """Given a list of complex numbers l and a complex number e, return the element e2 of l minimizing abs(e2-e)"""
        # find element in l that is closest to e for abs
        CC=ComplexField(500)
        r = l[0]
        for i in range(1,len(l)):
            if abs(CC(l[i]-e))<abs(CC(r-e)):
                r = l[i]
        return r

    @classmethod
    def select_closest_index(cls, l, e):
        """Given a list of complex numbers l and a complex number e, return the index i minimizing abs(l[i]-e)"""
        # find index of element in l that is closest to e for abs
        CC=ComplexField(500)
        r = 0
        for i in range(1,len(l)):
            if abs(CC(l[i]-e))<abs(CC(l[r]-e)):
                r = i
        return r

    @classmethod
    def is_clockwise(cls, l):
        """Given a list of complex numbers describing a convex polygon, return whether the points are clockwise."""
        CC=ComplexField(500)
        smally = min(l, key=lambda v:(CC(v).imag(), CC(v).real()))
        i = l.index(smally)
        n = l[i+1 if i+1<len(l) else 0]
        p = l[i-1]

        x1,x2,x3 = [v.real() for v in [p,smally,n]]
        y1,y2,y3 = [v.imag() for v in [p,smally,n]]

        M = matrix([[1,x1,y1],[1,x2,y2],[1,x3,y3]])
        if abs(CC(M.determinant()))<10e-7:
            logger.warning("cross product is very small, not certain about orientation")
        
        return CC(M.determinant())<0

    @classmethod
    def is_simple(cls, l):
        """Given a list of words l, return whether every word in the list consists of a single letter."""
        for w in l:
            if len(w.syllables()) != 1:
                return False
        return True
    
    @classmethod
    def letter(cls, w,i):
        """Given a word w and an integer i, yields the i-th letter of w."""
        return w.syllables()[i][0]**(w.syllables()[i][1]/abs(w.syllables()[i][1]))

    @classmethod
    def invert_morphism(cls, phi, ts = None):
        """Given an invertible free group morphism phi, computes its inverse. 
        Optionally, you can give generators ts of a subgroup (as a list of words) to compute the inverse of the restriction of phi on \\langle ts \\rangle (assuming it is invertible)."""
        # I have no idea if this terminates -- if it does it should not take very long
        pmax=20
        pcutoff=0
        if ts == None:
            ts = list(phi.domain().gens())
        while not Util.is_simple([phi(t) for t in ts]):
            managed = False
            for i, t in enumerate(ts):
                others = [t for j,t in enumerate(ts) if i != j]
                while len(phi(t).syllables())!=1:
                    options = [t*t2 for t2 in others if Util.letter(phi(t2),0) == Util.letter(phi(t),-1)**-1]
                    options += [t*t2**-1 for t2 in others if Util.letter(phi(t2**-1),0) == Util.letter(phi(t),-1)**-1]
                    options += [t2*t for t2 in others if Util.letter(phi(t2),-1) == Util.letter(phi(t),0)**-1]
                    options += [t2**-1*t for t2 in others if Util.letter(phi(t2**-1),-1) == Util.letter(phi(t),0)**-1]
                    if len(options)==0:
                        break
                    options = [o for o in options if phi(o)!=phi(1)]
                    options.sort(key = lambda w: len(phi(w).syllables()))
                    p=randint(1,pmax)
                    if len(phi(options[0]).syllables())< len(phi(t).syllables()):
                        t=options[0]
                        managed=True
                    else:
                        break
                ts[i] = t
            if not managed:
                for j, wi in enumerate(ts):
                    for i in range(j+1, len(ts)):
                        if Util.letter(phi(wi),0) != Util.letter(phi(wi),-1)**-1 and len(phi(wi).syllables())!=1:
                            while Util.letter(phi(wi),0) == Util.letter(phi(ts[i]),-1)**-1 or Util.letter(phi(wi),0) == Util.letter(phi(ts[i]),0):
                                if Util.letter(phi(wi),0) == Util.letter(phi(ts[i]),-1)**-1:
                                    ts[i] = ts[i]*wi
                                if Util.letter(phi(wi),0) == Util.letter(phi(ts[i]),0):
                                    ts[i] = wi**-1*ts[i]
                if pcutoff>5:
                    break
                pcutoff+=1
            else:
                pcutoff = max(0, pcutoff-1)
        tfin = [None]*len(ts)
        for t in ts:
            x, power = phi(t).syllables()[0]
            tfin[phi.codomain().gens().index(x)] = t**power
        return phi.codomain().hom(tfin)

    @classmethod
    def find_complement(cls, B):
        """Given an m x n integer valued matrix B with n>m, computes an (n-m) x n matrix A such that the matrix block_matrix([[A],[B]]) is invertible over the integers"""
        D, U, V = B.smith_form()
        quotient = identity_matrix(D.ncols())[D.nrows():]*V**-1
        assert block_matrix([[B],[quotient]]).det() in [-1,1], "cannot find complement, are you sure sublattice is primitive?"
        return quotient
    
    @classmethod
    def middle(self, w):
        """Given a word w of odd length 2n+1, yields the word consisting of the first n letters of w."""
        syls = w.syllables()
        syls = syls[:len(syls)//2]
        conj = w.parent(1)
        for l, p in syls:
            conj = conj*l**p
        return conj