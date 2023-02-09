# -*- coding: utf-8 -*-

import sage.all

from ore_algebra import *
from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.qqbar import AlgebraicField
from sage.rings.qqbar import QQbar
from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpfr import ComplexField
from sage.geometry.voronoi_diagram import VoronoiDiagram
from sage.graphs.graph import Graph
from sage.symbolic.constants import I
from sage.functions.other import arg
from delaunay_triangulation.triangulate import delaunay, Vertex
from sage.functions.other import ceil
from sage.functions.other import floor
from sage.symbolic.constants import pi
from sage.arith.misc import gcd
from sage.arith.misc import xgcd
from sage.plot.plot import list_plot
from sage.parallel.decorate import parallel
from sage.combinat.integer_vector import IntegerVectors
from sage.matrix.constructor import matrix


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
        return [ring.monomial(*m) for m in list(IntegerVectors(degree, ring.ngens()))]

    @classmethod
    def xgcd_list(cls, l):
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
        # find element in l that is closest to e for abs
        CC=ComplexField(500)
        r = l[0]
        for i in range(1,len(l)):
            if abs(CC(l[i]-e))<abs(CC(r-e)):
                r = l[i]
        return r

    @classmethod
    def select_closest_index(cls, l, e):
        # find index of element in l that is closest to e for abs
        CC=ComplexField(500)
        r = 0
        for i in range(1,len(l)):
            if abs(CC(l[i]-e))<abs(CC(l[r]-e)):
                r = i
        return r

    @classmethod
    def is_clockwise(cls, l):
        CC=ComplexField(500)
        # find index of element in l that is closest to e for abs
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

