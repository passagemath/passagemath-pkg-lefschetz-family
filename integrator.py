# -*- coding: utf-8 -*-

import sage.all

from numperiods import Family
from numperiods import Cohomology
from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.qqbar import AlgebraicField
from sage.rings.qqbar import QQbar
from sage.functions.other import ceil
from sage.functions.other import floor
from sage.functions.other import arg
from sage.functions.other import factorial
from sage.graphs.graph import Graph
from sage.symbolic.constants import I
from sage.rings.complex_double import CDF
from sage.matrix.constructor import matrix
from sage.arith.misc import xgcd
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import block_diagonal_matrix
from sage.matrix.special import zero_matrix
from sage.parallel.decorate import parallel
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from ore_algebra.analytic.monodromy import formal_monodromy
from ore_algebra.analytic.differential_operator import DifferentialOperator

from sage.misc.prandom import randint

from edges import Edges
from Util import Util
from Context import Context

import logging
import os
import time

logger = logging.getLogger(__name__)


class Integrator(object):
	"""This class takes a differential operator and a loop structure (either voronoi or delunay) and returns the numerical transition matrices along these paths
	"""
    def __init__(self, path_structure, operator):
        assert basepoint not in points

        self._operator =

        self.CC = ComplexField(500) # ultimately this should be dropped for certified precision





























