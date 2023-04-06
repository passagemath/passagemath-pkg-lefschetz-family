# -*- coding: utf-8 -*-

import sage.all

from ore_algebra import *

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.qqbar import QQbar
from sage.graphs.graph import Graph
from sage.symbolic.constants import I

from sage.rings.complex_mpfr import ComplexField
from sage.groups.free_group import FreeGroup
from sage.misc.flatten import flatten
from sage.schemes.curves.zariski_vankampen import followstrand


from Util import Util

import logging
from copy import copy

logger = logging.getLogger(__name__)


class RootsBraid(object):
    def __init__(self, P, edges, additional_points=[]):
        """P, a polynomial in two variables u and t.

        This class computes the braid group of roots (in t) of P(u) as u moves along a path
        """
        
        # assert P.is_homogeneous(), "nonhomogeneous defining polynomial"
        
        self.edges = edges
        self.P = P
        self.additional_points=additional_points
        self.npoints = self.P.degree(self.P.parent().gens()[1]) + len(self.additional_points)
        self._maximalstep = 1/1000

        self.freeGroup = FreeGroup(self.npoints)
        self.xs = list(self.freeGroup.gens())
        self.additional_points=additional_points


    @property
    def singularities(self):
        if not hasattr(self,'_singularities'):
            t = self.P.parent()('t')
            discrP = self.P.discriminant(t)
            Qu=PolynomialRing(QQ[I], 'u')
            self._singularities = Qu(discrP).roots(QQbar, multiplicities=False)
        return self._singularities
    

    def braid(self, e):
        if not hasattr(self,'_braid'):
            self._braid = [None]*len(self.edges)
            self._braidQ = [False]*len(self.edges)
        i=self.edge(e)
        if not self._braidQ[i]:
            CC=ComplexField(500)
            Qt=PolynomialRing(QQ[I], 't')
            u,t = self.P.parent()('u'),self.P.parent()('t')
            logger.info("Computing braid along edge %d."% i)
            roots=  Qt(self.P(u=e[0])).roots(QQbar, multiplicities=False)
            res=[]
            j=0
            for r in roots:   
                j+=1
                logger.info("Computing thread %d."% j)
                line = followstrand(self.P, [z.minpoly()(self.P.parent().gens()[1]) for z in self.additional_points], e[0], e[1],r, 50)
                res+=  [[[c[0], c[1]+I*c[2]] for c in line]]
            res.sort(key=lambda thread: (CC(thread[0][1].real()), CC(thread[0][1].imag())))
            self._braid[i] = res
            self._braidQ[i] = True

        return self._braid[i]

    def interpolate(self, path, t):
        if t==1:
            return path[-1][1]
        for i in range(len(path)-1,-1,-1):
            if t>=path[i][0]:
                break
        t0,x0, t1,x1 = path[i][0], path[i][1], path[i+1][0], path[i+1][1]
        return ((t1-t)*x0+ (t-t0)*x1)/(t1-t0)


    def minimal_cover_tree(self, section):
        CC=ComplexField(500)
        mtc=Graph(self.npoints+1) 
        edges = flatten([[(i,j) for i in range(j)] for j in range(self.npoints)], max_level=1) 
        edges.sort(key=(lambda e: Util.simple_rational(abs(CC(section[e[0]]-section[e[1]])), 10e-10))) # we sort edges by length
        for e in edges:
            if len(mtc.shortest_path(e[0], e[1]))==0:
                mtc.add_edge(e)
        # then we add the path to the basepoint
        vertices = [i for i in range(self.npoints)]
        vertices.sort(key=(lambda v: (section[v].real(), -section[v].imag()))) # the fixed basepoint is also specific to example
        mtc.add_edge([self.npoints,vertices[0]])
        return mtc

    def edge(self, e):
        return self.edges.index([e[0], e[1]])


    def transition_isomorphism(self,e1,e2):
        if not hasattr(self,'_transition_isomorphism'):
            self._transition_isomorphism = [[None for i in range(len(self.edges))] for j in range(len(self.edges))]
            self._transition_isomorphismQ = [[False for i in range(len(self.edges))] for j in range(len(self.edges))]
        
        i1 = self.edge(e1)
        i2 = self.edge(e2)
        if not self._transition_isomorphismQ[i1][i2]:
            logger.info("Computing transition isomorphism between edges %d and %d."% (i1,i2))

            braid1 = self.braid(e1)
            braid2 = self.braid(e2)

            section1 = self.braid_section(braid1, 1)
            section2 = self.braid_section(braid2, 0)

            perm = [Util.select_closest_index(section2,c) for c in section1]+[self.npoints] # this is fine because they are equal (although their presentation might differ)
            
            mtc1 = self.minimal_cover_tree(section1)
            mtc2 = self.minimal_cover_tree(section2)

            mtcn = Graph(self.npoints+1)
            for e in mtc1.edges():
                mtcn.add_edge((perm[e[0]], perm[e[1]]))

            oe1, oe2=self.ordered_edges(mtc1),self.ordered_edges(mtcn)
            perm_edge = [oe2.index(self.normalize_edge((perm[e[0]],perm[e[1]]))) for e in oe1]
            transition_iso = self.freeGroup.hom([self.xs[i] for i in perm_edge])
            self._transition_isomorphism[i1][i2] = self.braid_action(mtcn, mtc2, section2)*transition_iso
            self._transition_isomorphismQ[i1][i2] = True
        return self._transition_isomorphism[i1][i2]

    def braid_section(self, braid, t):
        section = []
        for thread in braid:
            section += [self.interpolate(thread, t)]
        return section+self.additional_points

    def isomorphisms(self, e):
        if not hasattr(self,'_isomorphisms'):
            self._isomorphisms=[None]*len(self.edges)
            self._isomorphismsQ=[False]*len(self.edges)
        i = self.edge(e)
        if not self._isomorphismsQ[i]:
            braid = self.braid(e)
            logger.info("Computing isomorphism of edge %d."% (i))
            ts=[]
            for thread in braid:
                for p in thread:
                    ts+=[p[0]]
            ts = list(set(ts))
            ts.sort()
            for t0,t1 in zip(ts[:-1], ts[1:]):
                while t1-t0>1.1*self._maximalstep:
                    t0+=self._maximalstep
                    ts.append(t0)
            ts.sort()
            sections = [self.braid_section(braid, t) for t in ts]

            mtcs = [self.minimal_cover_tree(section) for section in sections]

            iso = self.freeGroup.hom(self.xs)
            for k in range(len(mtcs)-1):
                mtc1 = mtcs[k]
                mtc2 = mtcs[k+1]
                if mtc1!=mtc2:
                    logger.info("Encountered distinct minimal covering trees between sections %d and %d (out of %d). Computing transition."% (k, k+1, len(ts)))
                    iso = self.braid_action(mtc1, mtc2, sections[k])*iso
                    iso = self.freeGroup.hom([iso(x) for x in self.xs])
            self._isomorphismsQ[i]=True
            self._isomorphisms[i]=iso
        return self._isomorphisms[i]
    
    def isomorphism_along_path(self,path):
        path_edges = [path[i:i+2] for i in range(len(path)-1)]
        iso = self.isomorphisms(path_edges[0])
        for i in range(len(path_edges)-1):
            iso = self.isomorphisms(path_edges[i+1])*self.transition_isomorphism(path_edges[i], path_edges[i+1])*iso
            iso = self.freeGroup.hom([iso(x) for x in self.xs])
        iso = self.transition_isomorphism(path_edges[-1], path_edges[0])*iso
        iso = self.freeGroup.hom([iso(x) for x in self.xs])
        return iso

    def normalize_edge(self, e):
        return (e[0], e[1]) if e[0]<=e[1] else (e[1], e[0])

    def braid_action(self,g1,g2, section):
        """ computes the isomorphism elements characterizing the change from graph g1 to graph g2
        """
        iso = self.freeGroup.hom(self.xs)
        removed_edges = []
        added_edges = []
        for e in g1.edges():
            if not g2.has_edge(e):
                removed_edges+=[e]
        for e in g2.edges():
            if not g1.has_edge(e):
                added_edges+=[e]
        # logger.info("Starting from graph with edges (%d, %d), (%d, %d), (%d, %d)."% tuple(flatten([[e[0], e[1]] for e in g1.edges()])))
        for e in removed_edges:
            # logger.info("Removing edge (%d,%d)."%(e[0], e[1]))
            # we add the edge to see what cycle appears
            gx = copy(g1)
            gx.delete_edge(e)

            if not e[0] in gx.connected_component_containing_vertex(self.npoints):
                e = (e[1], e[0])
            for ea in added_edges:
                ga = copy(gx)
                ga.add_edge(ea)
                if ga.connected_components_number()==1:
                    break
            ga = copy(g1)
            ga.add_edge(ea)
            assert ga.connected_components_number()==1
            # logger.info("Adding edge (%d,%d)."%(ea[0], ea[1]))
            ea = self.normalize_edge(ea)

            cycle = ga.cycle_basis()[0]
            # then we normalize the cycle (we want to start at the vertex from the edge that was deleted that's in the same component as the basepoint)

            ci = cycle.index(e[0])
            if cycle[ci+1 if ci+1<len(cycle) else 0] == e[1]:
                cycle.reverse()
            ci = cycle.index(e[0])
            cycle = cycle[ci:] + cycle[:ci]
            xmax=Util.simple_rational(max([s.real() for s in section]), 0.1)
            xmin=Util.simple_rational(min([s.real() for s in section]), 0.1)
            ymax=Util.simple_rational(max([s.imag() for s in section]), 0.1)
            clockwise = Util.is_clockwise([section[i] if i!=self.npoints else 2*xmin-xmax+ymax*I/5 for i in cycle])
            
            # logger.info("Clockwise cycle" if clockwise else "Counterclockwise cycle")

            beforeLoop = True
            beforeLink = True
            beforeEdges1 = []
            beforeEdges2 = []
            afterEdges = []
            link=min(cycle, key=lambda v:ga.distance(v, self.npoints))

            for j in range(len(cycle)-1):
                e2 = (cycle[j], cycle[j+1])
                if link==cycle[j]:
                    beforeLink = False
                if e2[0]>e2[1]:
                    e2=(e2[1],e2[0])
                if self.normalize_edge(e2) == ea:
                    beforeLoop = False
                    gx.add_edge(ea)
                else:
                    if beforeLink:
                        beforeEdges1+=[e2]
                    elif beforeLoop:
                        beforeEdges2+=[e2]
                    else:
                        afterEdges+=[e2]
            e = self.normalize_edge(e)


            res = []
            ini = self.ordered_edges(g1)
            fin = self.ordered_edges(gx)
            xr= self.xs[fin.index(ea)]

            for e2 in ini:
                if clockwise:
                    if e2 in beforeEdges1:
                        res+=[self.xs[fin.index(e2)]*xr]
                    elif e2 in beforeEdges2:
                        res+=[xr**-1*self.xs[fin.index(e2)]]
                    elif e2 in afterEdges:
                        res+=[self.xs[fin.index(e2)]**-1*xr]
                    elif e2 == e:
                        res+=[xr]
                    else:
                        res+=[self.xs[fin.index(e2)]]
                else:
                    if e2 in beforeEdges1:
                        res+=[xr*self.xs[fin.index(e2)]]
                    elif e2 in beforeEdges2:
                        res+=[self.xs[fin.index(e2)]*xr**-1]
                    elif e2 in afterEdges:
                        res+=[xr*self.xs[fin.index(e2)]**-1]
                    elif e2 == e:
                        res+=[xr]
                    else:
                        res+=[self.xs[fin.index(e2)]]

            iso = self.freeGroup.hom([self.freeGroup.hom(res)(iso(x)) for x in self.xs])
            g1=gx
        # logger.info("Ended with graph with edges (%d, %d), (%d, %d), (%d, %d)."% tuple(flatten([[e[0], e[1]] for e in g1.edges()])))
        assert g1==g2, "Did not recover correct graph"
        return iso

    def normalize_edge(self, e):
        return (e[0], e[1]) if e[0]<=e[1] else (e[1], e[0])

    def ordered_edges(self, gr):
        res = [(e[0], e[1]) for e in gr.edges()]
        for e in res:
            if e[0]>e[1]:
                e[0],e[1]=e[1],e[0]
        res.sort()
        return res






        