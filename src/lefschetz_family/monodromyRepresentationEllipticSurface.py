# -*- coding: utf-8 -*-

import sage.all

from .numperiods.familyNew import Family
from .numperiods.cohomology import Cohomology
from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar
from sage.matrix.constructor import matrix
from sage.arith.misc import xgcd
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.matrix.special import block_matrix
from sage.matrix.special import block_diagonal_matrix
from sage.matrix.special import zero_matrix
from sage.arith.functions import lcm
from sage.misc.misc_c import prod
from sage.misc.flatten import flatten

from sage.misc.prandom import randint

from .util import Util
from .context import Context
from .exceptionalDivisorComputer import ExceptionalDivisorComputer
from .monodromyRepresentation import MonodromyRepresentation
from .ellipticSingularity import EllipticSingularities


import logging
import time

logger = logging.getLogger(__name__)


class MonodromyRepresentationEllipticSurface(MonodromyRepresentation):

    @property
    def types(self):
        if not hasattr(self, "_types"):
            types = [EllipticSingularities.monodromy_class(M) for M in self.monodromy_matrices]
            types = [t + str(n) if t in ["I", "I*"] else t for t, _, n in types]
            self._types = types
        return self._types
    
    @property
    def monodromy_matrices_desingularisation(self):
        if not hasattr(self, '_monodromy_matrices_desingularisation'):
            I1_monodromy_matrices = []
            for M in self.monodromy_matrices:
                ty, base_change, nu = EllipticSingularities.monodromy_class(M)
                mats =  [base_change * M * base_change.inverse() for M in EllipticSingularities.fibre_confluence[ty][:-1]]
                mats += [base_change * EllipticSingularities.fibre_confluence[ty][-1] * base_change.inverse()] * nu
                mats = [M.change_ring(ZZ) for M in mats]
                Mtot = 1
                for M2 in mats:
                    Mtot = M2*Mtot
                assert Mtot == M
                I1_monodromy_matrices += [mats]

            self._monodromy_matrices_desingularisation = I1_monodromy_matrices

        return self._monodromy_matrices_desingularisation
    
    @property
    def self_intersection_section(self):
        if not hasattr(self, '_self_intersection_section'):
            chi = ZZ((len(self.extensions_desingularisation)+ZZ(4))/ZZ(12))
            self._self_intersection_section = -chi
        return self._self_intersection_section
    
    @property
    def add(self):
        if not hasattr(self, '_add'):
            self._add = 2
        return self._add
