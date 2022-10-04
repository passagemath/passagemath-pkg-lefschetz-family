# -*- coding: utf-8 -*-


from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpfr import ComplexField

class Context(object):

    def __init__(self,
            method=None,
            compute_periods=True,
            singular=False,
            debug=False,
            use_symmetry=True,
            nbits=200,
            depth=0
        ):
        r"""
        Lefschetz Family integration context

        Options:

        * ``method`` -- The way the paths are computed, either along a Voronoi diagram of the singularities ("voronoi"), or a Delaunay triangulation of the singularities ("delaunay"). Default is "voronoi"
        * ``compute_periods`` -- Whether the algorithm should compute periods of the variety, or stop at homology. Default is True.
        * ``singular`` -- Whether the input variety is expected to be singular. Default is False

        * (other options still to be documented...)
        """

        if not method in [None, "voronoi", "delaunay"]:
            raise ValueError("method", method)
        self.method = "voronoi" if method==None else method

        if not isinstance(compute_periods, bool):
            raise TypeError("compute_periods", type(compute_periods))
        self.compute_periods = compute_periods

        if not isinstance(singular, bool):
            raise TypeError("singular", type(singular))
        self.singular = singular
        
        if not isinstance(debug, bool):
            raise TypeError("debug", type(debug))
        self.debug = debug

        # if not isinstance(nbits, ): # what type is int ?
        #     raise TypeError("nbits", type(nbits))
        self.nbits = nbits

        if not isinstance(use_symmetry, bool):
            raise TypeError("use_symmetry", type(use_symmetry))
        self.use_symmetry = use_symmetry

        # if not isinstance(depth, int):
        #     raise TypeError("depth", type(depth))
        # self.depth = depth

        self.CBF = ComplexBallField(500)
        self.CF = ComplexField(500)
        self.depth = depth

dctx = Context() # default context
