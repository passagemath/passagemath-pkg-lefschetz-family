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
from sage.plot.plot import list_plot
from sage.parallel.decorate import parallel
from sage.combinat.integer_vector import IntegerVectors

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
