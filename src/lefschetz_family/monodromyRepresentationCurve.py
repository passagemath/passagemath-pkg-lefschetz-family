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

from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ

from .monodromyRepresentation import MonodromyRepresentation

import logging
import time

logger = logging.getLogger(__name__)


class MonodromyRepresentationCurve(MonodromyRepresentation):

    def desingularise_matrix(self, M):
        if M==1:
            return []
        if (M-1).rank() != 1:
            raise Exception("Unknown singular fibre type")
        v = (M-1).image().gen(0)
        n = gcd(v)
        decomposition = [(M-1)/n + 1] * n
        decomposition = [M.change_ring(ZZ) for M in decomposition]
        return decomposition
    
    @property
    def add(self):
        if not hasattr(self, '_add'):
            self._add = 0
        return self._add