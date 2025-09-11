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

try:
    import sage.all
except ImportError:
    import sage.all__sagemath_modules

from ore_algebra import *

from sage.modules.free_module_element import vector
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.matrix.special import identity_matrix

from sage.functions.other import floor
from sage.functions.other import ceil
from sage.arith.misc import gcd

from sage.symbolic.relation import solve
from sage.symbolic.ring import SR

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
    
    def desingularise_matrix(self, M):
        ty, base_change, nu = EllipticSingularities.monodromy_class(M)
        decomposition = EllipticSingularities.fibre_confluence[ty]
        mats =  [base_change * M * base_change.inverse() for M in decomposition[:-1]]
        mats += [base_change * decomposition[-1] * base_change.inverse()] * nu
        mats = [M.change_ring(ZZ) for M in mats]
        return mats

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