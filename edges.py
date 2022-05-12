# -*- coding: utf-8 -*-

import sage.all

from sage.rings.qqbar import AlgebraicField
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpfr import ComplexField
from sage.graphs.graph import Graph
from sage.symbolic.constants import I
from sage.functions.other import arg
from delaunay_triangulation.triangulate import delaunay, Vertex
from sage.functions.other import floor
from sage.functions.other import ceil

from sage.geometry.voronoi_diagram import VoronoiDiagram



from Util import Util

#   TODO : overhaul of the Delaunay triangulation algorithm so that vertices of the graph are not algebraic/rational numbers but integer, using an index system
#
#

class Edges(object):

    CC=ComplexField(500) # choice is arbitrary, need a better way

    @classmethod
    def minimal_cover_tree(cls, graph, points, req_points=[]):
        """ Given a planar graph in the complex plane, computes a subtree of this graph, favoritizing smallest edges
        Optionnally, if a list of points is given, computes the smallest subtree containing all the points in the list.

        """
        CC=cls.CC

        sorted_edges = list(graph.edges())
        sorted_edges.sort(key=lambda e: abs(CC(points[e[0]]-points[e[1]])))

        G=Graph() # there should be a better way
        G2=Graph()
        for e in sorted_edges:
            G2.add_edge(e)
            if len(G2.cycle_basis())==0:
                G.add_edge(e)
            G2 = G.copy()
            if len(G.vertices()) == len(graph.vertices()):
                break

        if len(req_points)==0 or len(req_points)==1:
            return G

        cover_tree=Graph()
        

        cover_tree.add_edge(G.edges(vertices=[req_points[0]], labels=False)[0])
        for p in req_points[1:]:
            assert G.has_vertex(p), "Graph does not contain this point"
            sp=G.shortest_path(req_points[0], p)
            for i in range(len(sp)-1):
                cover_tree.add_edge(sp[i], sp[i+1])
        cover_tree.remove_multiple_edges()

        return cover_tree

    @classmethod
    def angle_positive(cls, a,b,c):
        CC= cls.CC
        r1 = (c-b)
        u=r1.real()
        v=r1.imag()
        r2 = (a-b)
        x = r2.real()
        y = r2.imag()
        im = CC(v*x-u*y)
        return im<0 

    @classmethod
    def angle_zero(cls, a,b,c):
        CC= cls.CC
        r1 = (c-b)
        u=r1.real()
        v=r1.imag()
        r2 = (a-b)
        x = r2.real()
        y = r2.imag()
        real = CC(u*x+v*y)
        im = v*x-u*y
        return cls.is_zero(im) and real>0 

    @classmethod
    def break_edges(cls, graph, edge, values):
        """ Gives a path in graph between both vertices of edge, taking into account the arrival angle
        """
        CC = cls.CC

        sp = graph.shortest_path(edge[0], edge[1], by_weight=True, weight_function=(lambda e: abs(CC(values[e[0]])-CC(values[e[1]]))))

        a,b,c,d = values[sp[-2]], values[edge[1]], values[edge[0]], values[sp[1]]

        if ((not cls.angle_positive(a, b, b-1)) and cls.angle_positive(c, b, b-1) and cls.angle_positive(c, b, a)):
            sp = sp + ["Rd"]
        elif (cls.angle_positive(a, b, b-1) and cls.angle_positive(b-1, b, c) and cls.angle_positive(a, b, c)):
            sp = sp + ["Ri"]
        elif (cls.angle_zero(a, b, b-1) and cls.angle_positive(b-1, b, c)):
            sp = sp + ["Ri"]
            
        if ((not cls.angle_positive(d, c, c-1)) and cls.angle_positive(b, c, c-1) and cls.angle_positive(b, c, d)):
            sp = sp + ["Ri"]
        elif (cls.angle_positive(d, c, c-1) and cls.angle_positive(c-1, c, b) and cls.angle_positive(d, c, b)):
            sp = sp + ["Rd"]
        elif (cls.angle_zero(d, c, c-1) and cls.angle_positive(c-1, c, b)):
            sp = sp + ["Rd"]

        return sp

    @classmethod
    def reverse(cls, path):
        """ reverses a path (taking into account Ri and Rd)
        """
        res = []
        for p in path:
            if p=="Ri":
                res = ["Ri"] + res
            elif p=="Ri":
                res = ["Rd"] + res
            else:
                res = [p] + res
        return res
    
    @classmethod
    def angle(cls, a, b, c):
        CC = cls.CC
        return arg(-(CC(a)-CC(b))/(CC(c)-CC(b)))


    @classmethod
    def sort_neighbors_delaunay(cls, p, neighbors, alg_points):
        res = list(neighbors)
        res.sort(key = lambda x: -cls.angle(alg_points[p]-1, alg_points[p], alg_points[x]))
        if Edges.angle_zero(alg_points[p]-1, alg_points[p], alg_points[res[0]]):
            res = res[1:]+[res[0]]
        return res

    @classmethod
    def _circuit_delaunay_rec(cls, neighbors, basepoint, points, actual, previous=None, start=False):
        if start:
            return cls._circuit_delaunay_rec(neighbors, basepoint, points, neighbors[actual][0], actual)
        assert previous!=None, "circuit badly intialised"

        if actual==basepoint:
            return []
        neigh = neighbors[actual]
        i = neigh.index(previous)
        if i == len(neigh)-1:
            if actual in points:
                return [points.index(actual)] + cls._circuit_delaunay_rec(neighbors, basepoint, points, neigh[0], actual)
            return cls._circuit_delaunay_rec(neighbors, basepoint, points, neigh[0], actual)
        return cls._circuit_delaunay_rec(neighbors, basepoint, points, neigh[i+1], actual)

    @classmethod
    def circuit_delaunay(cls, T, points, basepoint, alg_points):
        """ Returns a pair paths, order where
        - paths is a list of paths ps in the graph from basepoint to each element of points 
        - order is a permutation of [1, ..., len(points)] such that the composition ps[order[len(points)-1] * ... * ps[order[0]] is the simple trigonometric loop around all points of points
        
        Note that cuts in ore_algebra are always along the negative real line

        """

        neighbors={
            "blank":None
        }
        for v in T.vertices():
            neighbors[v] = cls.sort_neighbors_delaunay(v, T.neighbors(v), alg_points)
        order=cls._circuit_delaunay_rec(neighbors, basepoint, points, basepoint, start=True)
        
        
        paths=[T.shortest_path(basepoint, p) for p in points]

        return paths, order

    @classmethod
    def is_zero(cls, x, eps=2**-50):
        return abs(cls.CC(x))<eps

    @classmethod        
    def delaunay(cls, points):
        """Computes the Delaunay triangulation of a cloud of points in the complex plane

        """
        alg_points = [QQbar(p) for p in points]

        CC = cls.CC
        r  = min([abs(CC(x-y)) for x in alg_points for y in alg_points if x!=y])

        Qpoints = [Util.simple_rational(p.real(), r/3) + I*Util.simple_rational(p.imag(), r/3) for p in alg_points]
        dictionary = {Qpoints[i]:i for i in range(len(alg_points))}

        # we deal with rational points because equality is faster
        triangles = delaunay(vertices=[Vertex(x=Qp.real(), y=Qp.imag()) for Qp in Qpoints])
        G=Graph()

        for t in triangles:
            for e in t.edges:
                e1=e.a.to_tuple[0]+e.a.to_tuple[1]*I
                e2=e.b.to_tuple[0]+e.b.to_tuple[1]*I
                G.add_edge((dictionary[e1],dictionary[e2]))

        return  G

    @classmethod
    def closest(cls, p, l): # not used
        """ Finds the element closest to p in l for the euclidean norm in the complex plane
        """
        res=l[0]
        for e in l:
            if abs(p-e) < abs(p-res):
                res = e
        return res
    
    @classmethod
    def voronoi_loops(cls, points, basepoint):

        vd = cls.voronoi(points, basepoint)

        gr, loops = cls.loops(points, vd, basepoint)

        cover_tree = cls.subtree(basepoint, [l[0] for l in loops], gr)

        sps=[cover_tree.shortest_path(basepoint, l[0]) for l in loops]

        order = cls.voronoi_circuit(points, basepoint, cover_tree, loops)

        return sps, loops, order

    @classmethod
    def voronoi(cls, points, basepoint, shift=1):

        CC = cls.CC

        r = min([min([ abs(CC(points[i]-points[j])) for j in range(i+1, len(points))]) for i in range(len(points)-1)])/10

        maxx=ceil(max([p.real() for p in points])+shift)
        minx=floor(min([p.real() for p in points])-shift)
        maxy=ceil(max([p.imag() for p in points])+shift)
        miny=floor(min([p.imag() for p in points])-shift)
               
        rootapprox = set((QQ(Util.simple_rational(p.real(), r)), QQ(Util.simple_rational(p.imag(), r))) for p in points)
        rootapprox.add((QQ(basepoint.real()), QQ(basepoint.imag())))
        rootapprox.add((QQ(maxx), QQ(0)))
        rootapprox.add((QQ(0), QQ(maxy)))
        rootapprox.add((QQ(0), QQ(miny))) 
        rootapprox.add((QQ(maxx), QQ(miny))) 
        rootapprox.add((QQ(maxx), QQ(maxx))) 
        rootapprox.add((QQ(minx), QQ(minx))) 
        rootapprox.add((QQ(minx), QQ(maxx))) 

        return VoronoiDiagram(rootapprox)


    @classmethod
    def loops(cls, points, vd, basepoint):

        CC= cls.CC
        r = min([min([ abs(CC(points[i]-points[j])) for j in range(i+1, len(points))]) for i in range(len(points)-1)])/10
        
        gr = Graph()
        for pt_, reg in vd.regions().items():
            pt = tuple(pt_)
            pt = pt[0] + I*pt[1]
            for edge in reg.bounded_edges():
                u = Util.simple_rational(edge[0][0],r)+I*Util.simple_rational(edge[0][1],r)
                v = Util.simple_rational(edge[1][0],r)+I*Util.simple_rational(edge[1][1],r)
                if u!=v:
                    gr.add_edge(u, v)

            if pt == basepoint:
                for v in reg.vertices():
                    u = Util.simple_rational(v[0],r)+I*Util.simple_rational(v[1],r)
                    gr.add_edge(pt, u, 1.0)
                    break
        gr.remove_multiple_edges()


        regions = [] # lists all regions containing critical points
        for p in [[CC(s).real().simplest_rational(), CC(s).imag().simplest_rational()] for s in points]:
            for pt_, reg in vd.regions().items():
                if reg.contains(p):
                    regions.append(reg)
                    break

        assert len(regions) == len(points)

        loops=[]
        for reg in regions:
            loops+=[[]]
            edges = [[Util.simple_rational(e[0][0],r)+Util.simple_rational(e[0][1],r)*I, Util.simple_rational(e[1][0],r)+Util.simple_rational(e[1][1],r)*I] for e in reg.bounded_edges()]
            loops[-1] += edges[0]
            edges.pop(0)
            i=0
            while len(edges)>1:
                if edges[i][0]==loops[-1][-1]:
                    loops[-1]+=[edges[i][1]]
                    edges.pop(i)
                    i=0
                elif edges[i][1]==loops[-1][-1]:
                    loops[-1]+=[edges[i][0]]
                    edges.pop(i)
                    i=0
                else:
                    i+=1

        # remove consecutive doubles
        for j in range(len(loops)):
            loop = loops[j]
            i=0
            while i<len(loop)-1:
                if loop[i]==loop[i+1]:
                    loop = loop[:i+1]+loop[i+2:]
                else:
                    i+=1
            loops[j]=loop



        # this ensures the loops are trigonometric-wise around the critical points
        for loop in loops:
            y = loop[0]
            i=0
            for k in range(len(loop)):
                x=loop[k]
                if x.imag()>y.imag() or (x.imag()==y.imag() and x.real()>y.real()):
                    y=x
                    i=k
            if i==0:
                A=loop[0]
                B=loop[-1]
                C=loop[1]
            elif i==len(loop)-1:
                A=loop[-1]
                B=loop[-2]
                C=loop[0]
            else:
                A=loop[i]
                B=loop[i-1]
                C=loop[i+1]
            if not((B.real()-A.real())*(C.imag()-A.imag())-(B.imag()-A.imag())*(C.real()-A.real())<0):
                loop = loop.reverse()

        return gr, loops    
    
    @classmethod
    def subtree(cls, basepoint, l, gr):
        cover_tree=Graph()
        cover_tree.add_edge(gr.edges(basepoint)[0])
        gr2=gr.copy()
        for v in l:
            if gr2.has_vertex(v):
                minDist = gr2.shortest_path_length(basepoint, v)
                for e in cover_tree.vertices():
                    sp=gr.shortest_path(e, v)
                    if len(sp)==minDist:
                        break
                for p in sp:
                    if gr2.has_vertex(p):
                        gr2.merge_vertices([basepoint, p])
                for i in range(len(sp)-1):
                    cover_tree.add_edge(sp[i], sp[i+1])
        cover_tree.remove_multiple_edges()

        # making sure it's a tree
        basis = cover_tree.cycle_basis()
        while len(basis) !=0:
            cover_tree.delete_edge(basis[0][0], basis[0][1])
            basis = cover_tree.cycle_basis()

        sps=[cover_tree.shortest_path(basepoint, v) for v in l]
        return cover_tree


    @classmethod
    def voronoi_circuit_rec(cls, neighbors, basepoint, loops, actual, previous=None, start=False):
        if start:
            return cls.voronoi_circuit_rec(neighbors, basepoint, loops, neighbors[actual][0], actual)
        if actual == basepoint:
            return []
        assert previous!=None, "badly intialised"

        neigh = neighbors[actual]
        i = neigh.index(previous)
        if i == len(neigh)-1:
            n = neigh[0]
        else:
            n = neigh[i+1]
        order = []
        for i in range(len(loops)):
            if loops[i][0] == actual and (previous == n or cls.angle(previous, actual, n) < cls.angle(previous, actual, loops[i][1]) or cls.angle_zero(previous, actual, loops[i][1])):
                order+=[i]
        order.sort(key = lambda i: -cls.angle(previous, actual,loops[i][1]))

        
        return order+ cls.voronoi_circuit_rec(neighbors, basepoint, loops, n, actual)


    @classmethod
    def sort_neighbors_voronoi(cls, p, neighbors):
        res = list(neighbors)
        res.sort(key = lambda x: -cls.angle(p-1, p, x))
        if Edges.angle_zero(p-1, p, res[0]):
            res = res[1:]+[res[0]]
        return res    


    @classmethod
    def voronoi_circuit(cls, points, basepoint, T, loops):

        neighbors={
            "blank":None
        }
        for v in T.vertices():
            neighbors[v] = cls.sort_neighbors_voronoi(v, T.neighbors(v))
        order=cls.voronoi_circuit_rec(neighbors, basepoint, loops, basepoint, start=True) 

        return order



