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
from sage.matrix.special import zero_matrix
from sage.parallel.decorate import parallel
from sage.arith.misc import gcd
from ore_algebra.analytic.monodromy import formal_monodromy

from edges import Edges
from Util import Util
from Context import Context

import logging
import os

logger = logging.getLogger(__name__)


class LefschetzFamily(object):
    def __init__(self, P, **kwds):
        """P, a homogeneous polynomial defining a smooth hypersurface X in P^{n+1}.

        This class aims at computing an effective basis of the homology group H_n(X), 
        given as lifts of paths through a Lefschetz fibration.
        """

        logger.info("Computing homology of hypersurface defined by %s"% str(P))
        
        self.ctx = Context(**kwds)
        
        assert P.is_homogeneous(), "nonhomogeneous defining polynomial"
        
        self.degree = P.degree()
        self.P = P
        self._nvars = P.parent().ngens()
        self.dim = self._nvars-2
        self._R = P.parent()
        self._vars = [v for v in self._R.gens()]
        
        
        logger.info("Computing homology of variety of dimension %d." % self.dim)
        
        if self.dim==0:
            self.cohomology = Cohomology(self.P)
            affineR = PolynomialRing(QQbar, 'X')
            affineProjection= self._R.hom([affineR.gens()[0],1], affineR)
            self.homology = [e[0] for e in affineProjection(self.P).roots()]
            
            self.period_matrix = matrix([self._residue_form(affineProjection(b), affineProjection(self.P), (b.degree()+len(self._R.gens()))//P.degree(), self.homology) for b in self.cohomology.basis()]).change_ring(self.ctx.CBF)
            self.intersection_product=identity_matrix(self.degree)
        
        else:
            self.cohomology = Cohomology(self.P)
            self.fibration = self._compute_fibration()
            self.critical_points = self._compute_critical_points()
        
            self._S = PolynomialRing(PolynomialRing(QQ, self._vars[:-1]), 't')
            t=self._S.gens()[0]
            form = self.fibration[1][:-1].dot_product(vector(self._S.base_ring().gens())) + t*self.fibration[0][:-1].dot_product(vector(self._S.base_ring().gens()))
            denom = self.fibration[1][-1]*t+self.fibration[0][-1]
            self._RtoS = self._R.hom([denom**self.degree*v for v in self._S.base_ring().gens()]+[t*form], self._S)
            self.family = Family(self._RtoS(self.P))


            # thimbles are given as a list of (p, l) where p is an element of H_{n-1}(X_b) and l a path 
            # in P^1 \ Sigma such that tau_l(p) = Delta is the thimble
            self.thimbles = []
            self.monodromy_matrices =[]
            self.vanishing_cycles =[]
            self.perm_cycles = []
            list_thimbles = self._compute_thimbles() # the basepoint is defined here
            for (d, l, M) in list_thimbles:
                self.thimbles+=[(d, l)]
                self.perm_cycles += [d]
                if not self.ctx.debug:
                    self.vanishing_cycles+=[(M-1)*d / gcd((M-1)*d)]
                self.monodromy_matrices += [M]
            
            # homology cycles are given as a vector corresponding to coordinates of thimbles
            if not self.ctx.debug:
                self.homology = self._compute_homology()
                self.intersection_product=self._compute_intersection_product()
                if self.ctx.compute_periods:
                    self.period_matrix = self._compute_periods()
        logger.info("Computation of homology of variety of dimension %d completed." % self.dim)
    
    def _compute_fibration(self):
        # TODO: find good fibration to separate critical points. Ideally, minimize standard deviation of distance between singularities (?)
        # best_fibration=None
        # best_distance=0 
        # for i in range(100):
        #     self.fibration = (vector([0]*(self.dim+1)+[1]), vector([ZZ.random_element(-100,100) for i in range(self.dim+1)]+[0]))
        #     self.critical_points = self._compute_critical_points()
        #     distance = self._compute_distance(self.critical_points)
        #     if best_fibration==None or best_distance<distance:
        #         best_fibration = self.fibration
        #         best_distance = distance
        # self.fibration = best_fibration

        if self.dim==1:
            return (vector([0,0,1]), vector([2,5,0]))
        if self.dim==2:
            return (vector([0,0,0,1]),vector([2,5,7,0]))
        if self.dim==3:
            return (vector([0,0,0,0,1]),vector([2,5,7,11,0]))

    # def _compute_distance(self, points):
    #     return min([min([abs(ComplexField(50)(p1-p2)) for p2 in points if p2!=p1]) for p1 in points])
        

    def _compute_critical_points(self):
        """Returns a list of the critical points of the projection map, as algebraic complex numbers
        """
        forms=[v.dot_product(vector(self._vars)) for v in self.fibration]
        f=forms[0]/forms[1]
        S = PolynomialRing(QQ, self._vars+['k','t'])
        k,t= S.gens()[-2:]
        eqs = [self.P, self._vars[-1]-1,t*forms[1]-forms[0]]+[(f.derivative(var).numerator()-k*self.P.derivative(var)*f.derivative(var).denominator()) for var in self._vars]

        ideal = S.ideal(eqs).elimination_ideal(S.gens()[:-1])
        Qt = PolynomialRing(QQ, 't')
        return [e[0] for e in Qt(ideal.groebner_basis()[0]).roots(AlgebraicField())]

               
    def _compute_thimbles(self):
        """ Returns a list of triples (v, path, M) for each critical point x in self.critical_points, where 
            -- `path` is a loop around x, 
            -- `M` is the matrix of the action of monodromy along path on the homology of self.fiber
            -- `v` is a vector such that Mv-v is the vanishing cycle of x
        """
        
        logger.info("Computing thimbles of dimension %d ..." % self.dim)
               
        omega = vector([0]*(len(self.family.basis)-1)+[1]) # choice of cyclic vector
        L = self.family.picard_fuchs_equation(omega)
        assert L.order()== len(self.family.basis)

        paths = self._compute_paths()

        transition_matrices= self.integrate(L)

        logger.info("Computing periods of fiber.")
        self.evaluate_at_basepoint = self._S.hom([self.basepoint], self._S.base_ring())
        self.fiber = LefschetzFamily(self.evaluate_at_basepoint(self._RtoS(self.P)), method=self.ctx.method, nbits=self.ctx.nbits, depth=self.ctx.depth+1)


        n = len(self.fiber.homology)
        r = len(self.critical_points)   

        if self.dim%2==1:
            invariant_vector=vector([1]*n)
            proj = block_matrix([[identity_matrix(n-1), matrix([[-1]]*(n-1))]],subdivide=False)
            lift = block_matrix([[identity_matrix(n-1)], [matrix([[0]*(n-1)])]],subdivide=False)
            basis_change = block_matrix([[lift, matrix([-invariant_vector]).transpose()]],subdivide=False)
               
        derivatives = [omega*vector(self.family.basis)]
        for k in range(n-1 if self.dim%2==0 else n-2):
            derivatives += [self._derivative(derivatives[-1], self._RtoS(self.P))/(k+1)] 
        derivatives_coordinates = matrix([self.fiber.cohomology.coordinates(self.evaluate_at_basepoint(d)) for d in derivatives])
        
        if self.dim%2==1:
            initial_conditions = derivatives_coordinates*self.fiber.period_matrix*block_matrix([[identity_matrix(n-1)], [matrix([[0]*(n-1)])]])
        else:
            initial_conditions = derivatives_coordinates*self.fiber.period_matrix
        
        logger.info("Computing monodromy matrices")

        Mtot=1 # TODO: delete these lines
        for M in transition_matrices:
            Mtot=M*Mtot

        Ms = [(initial_conditions**(-1)*M*initial_conditions) for M in transition_matrices]
        if not self.ctx.debug:
            Ms = [M.change_ring(ZZ) for M in Ms]

        for i in range(r):
            M=Ms[i]
            if not self.ctx.singular and not self.ctx.debug:
                assert (M-1).rank()==1, "If M is a monodromy matrix around a single critical point, M-I should have rank 1"
            if self.dim%2==1:
                prod= invariant_vector*self.fiber.intersection_product*(lift*M*proj-identity_matrix(n))*lift/n
                Ms[i] = basis_change.inverse()*block_matrix([[M,zero_matrix(n-1,1)],[matrix([prod]),matrix([[1]])]])*basis_change

        vs=[]
        for M in Ms:
            for i in range(n):
                v = vector([0 if j!=i else 1 for j in range(n)]) # I am not sure one of these necessarily maps to a generator of the image
                if (M-1)*v!=0:
                    vs+=[v]
                    break
        
        logger.info("Found %d thimbles in dimension %d." % (r, self.dim))

        if not self.ctx.debug:
            Mtot = 1
            for M in Ms:
                Mtot = M*Mtot
            assert Mtot == identity_matrix(len(self.fiber.homology)), "Monodromy around infinity is nontrivial, most likely due to a mistake while computing a basis of homotopy."
               
        return [(vs[i], paths[i], Ms[i]) for i in range(r)]
               
    def _compute_homology(self):
        
        logger.info("Gluing thimbles together in dimension %d ..." % self.dim)
        
        delta = matrix(self.vanishing_cycles)
        kerdelta= delta.kernel()

        r = len(self.monodromy_matrices)
        n = len(self.fiber.cohomology.basis())
        
        M=1
        phi=[]
        for i in range(r):
            tempM=(self.monodromy_matrices[i]-1)*M
            phi+=[[c/self.vanishing_cycles[i] for c in tempM.columns()]]
            M=self.monodromy_matrices[i]*M
        phi = matrix(phi).transpose().change_ring(ZZ)
        assert M == identity_matrix(len(self.fiber.homology)), "monodromy around infinity is not trivial, either there is a critical point at infinity or the paths are not ordered properly"
               
        imphi = phi.image()
               
        B=matrix(kerdelta.gens()).solve_left(matrix(imphi.gens()))
        Brows=B.row_space()
               
        compl = [[0 for i in range(Brows.degree())]]
        rank=Brows.dimension()
        N=0
        for i in range(Brows.degree()):
            v=[1 if j==i else 0 for j in range(Brows.degree())]
            M=block_matrix([[B],[matrix(compl)],[matrix([v])]],subdivide=False)
            if rank+N+1==M.rank():
                compl += [v]
                N+=1
            if rank+N == Brows.degree():
                break
        quotient_basis=matrix(compl[1:])
        result = (quotient_basis*matrix(kerdelta.gens())).rows()
        logger.info("Done. Homology has rank %d in dimension %d." % (len(result), self.dim))
               
        return result
    
    # Integration methods

    def integrate(self, L):
        if self.ctx.method == "delaunay":
            sps, edges, values = self._break_edges(L)
            logger.info("Computing numerical transition matrices along edges (%d edges total)."% len(edges))
            ntms = [[r[0][0][1], r[1]] for r in list(LefschetzFamily._compute_transition_matrix_delaunay([(L,[e[0], e[1]], values, self.ctx.nbits) for e in edges]))]
            transition_matrices = []
            for j in range(len(self.critical_points)):
                M = self._reconstruct_path_delaunay(sps[j], ntms, L, values).change_ring(self.ctx.CBF)
                transition_matrices+=[M**-1*formal_monodromy(L, self.critical_points[j], ring=self.ctx.CBF)*M]
        elif self.ctx.method == "voronoi":
            N = len(self.edges)
            logger.info("Computing numerical transition matrices along edges (%d edges total)."% N)
            ntms = [[r[0][0][2], r[1]] for r in list(LefschetzFamily._compute_transition_matrix_voronoi([([i,N],L,[e[0], e[1]], self.ctx.nbits) for i, e in list(enumerate(self.edges))]))]
            transition_matrices = []
            for j in range(len(self.critical_points)):
                path =  self._sps[j] + self._loops[j] +list(reversed(self._sps[j]))
                M = self._reconstruct_path_voronoi(Util.simplify_path(path), ntms, L)
                transition_matrices+=[M]
        return transition_matrices

    def _reconstruct_path_delaunay(self, sp, ntms, L, values):
        if len(sp)<=1:
            return 1
        e = sp[:2]
        if e[1] == "Rd":
        	return self._reconstruct_path_delaunay([sp[0]]+sp[2:], ntms, L, values)*formal_monodromy(L, values[e[0]], ring=self.ctx.CBF)
        if e[1] == "Ri":
        	return self._reconstruct_path_delaunay([sp[0]]+sp[2:], ntms, L, values)*formal_monodromy(L, values[e[0]], ring=self.ctx.CBF)**-1
        if e[0] == "Rd":
        	return self._reconstruct_path_delaunay(sp[1:], ntms, L, values)*formal_monodromy(L, values[e[1]], ring=self.ctx.CBF)
        if e[0] == "Ri":
        	return self._reconstruct_path_delaunay(sp[1:], ntms, L, values)*formal_monodromy(L, values[e[1]], ring=self.ctx.CBF)**-1
        if e[0]==e[1]:
            return self._reconstruct_path_delaunay(sp[1:], ntms, L, values)
        for e2,M in ntms:
            if e2==e:
                return self._reconstruct_path_delaunay(sp[1:], ntms, L, values)*M
            if list(reversed(e2))==e:
                return self._reconstruct_path_delaunay(sp[1:], ntms, L, values)*M**-1
        print(e)
        raise Exception("unknown edge") 

    @classmethod
    def _reconstruct_path_voronoi(self, sp, ntms, L):
        if len(sp)<=1:
            return 1
        e = sp[:2]
        if e[0]==e[1]:
            return self._reconstruct_path_voronoi(sp[1:], ntms, L)
        for e2,M in ntms:
            if e2==e:
                return self._reconstruct_path_voronoi(sp[1:], ntms, L)*M
            if list(reversed(e2))==e:
                return self._reconstruct_path_voronoi(sp[1:], ntms, L)*M**-1
        print(e)
        raise Exception("unknown edge") 
               
    def _break_edges(self, L):
        logger.info("Adapting integration paths to delaunay triangulation of singularities of operator")

        sings = [self.basepoint] + self.critical_points + [x for x in L.leading_coefficient().roots(QQbar, multiplicities=False) if not x in self.critical_points]
        Delaunay = Edges.delaunay(sings)

        logger.info("Breaking edges (%d total)"% len(self.edges))
        paths=[]
        for e in self.edges:
            paths+=[Edges.break_edges(Delaunay, e, sings)]

        edges = []
        for p in paths:
            for i in range(len(p)-1):
                if p[i]!= "Ri" and p[i]!= "Rd" and p[i+1]!= "Ri" and p[i+1]!= "Rd":
                    if not [p[i], p[i+1]] in edges and not [p[i+1], p[i]] in edges:
                        edges += [[p[i], p[i+1]]]
        
        logger.info("Reconstructing paths with broken edges")
        sps = []    
        for sp in self._sps:
            sps+=[[]]
            for i in range(len(sp)-1):
                e = sp[i:i+2]
                for i in range(len(self.edges)):
                    e2 = self.edges[i]
                    if e[0]==e2[0] and e[1] == e2[1]:
                        sps[-1]+=paths[i]
                        break
                    if e[0]==e2[1] and e[1] == e2[0]:
                        sps[-1]+=Edges.reverse(paths[i])
                        break

        return sps, edges, sings
        


    def _compute_periods(self):
        
        logger.info("Computing periods in dimension %d." % self.dim)
        integrated_thimbles=[]
        
        N=len(self.cohomology.basis())
        n=len(self.fiber.homology)
        r=len(self.thimbles)
        
        Ls = self._compute_picard_fuchs()
        # self._Ls = [L*L.parent().gens()[0] for L in Ls]

        for i in range(N):
            logger.info("Computing integral of form along thimbles [%d/%d]"% (i+1, N))
            L = Ls[i]*Ls[i].parent().gens()[0]
            w = self.cohomology.basis()[i]
            wt = self._restrict_form(w)
            derivatives = [self._RtoS(0), wt]
            for k in range(n-1 if self.dim%2==0 else n-2):
                derivatives += [self._derivative(derivatives[-1], self._RtoS(self.P))] 
            derivatives_coordinates, denom = self.family.coordinates(derivatives) # this is not optimal - it would be better to compute all coordinates together
            
            integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(n+1 if self.dim%2==0 else n)])
            initial_conditions = integration_correction* derivatives_coordinates(self.basepoint)/denom(self.basepoint)*self.fiber.period_matrix
            
            transition_matrices = self.integrate(L)

            integrated_thimbles+=[[(transition_matrices[j]*initial_conditions*self.perm_cycles[j])[0] for j in range(r)]]
        
        logger.info("Integration of forms along thimbles computed in dimension %d." % self.dim)
        
        self._integrated_thimbles = integrated_thimbles
        
        return matrix(integrated_thimbles)*matrix(self.homology).transpose()
    
    def _compute_intersection_product(self):
        r=len(self.thimbles)
        inter_prod_thimbles = matrix([[self._compute_intersection_product_thimbles(i,j) for j in range(r)] for i in range(r)])
        return (matrix(self.homology)*inter_prod_thimbles*matrix(self.homology).transpose()/2).change_ring(ZZ)
        
    def _compute_intersection_product_thimbles(self,i,j):
        vi = self.thimbles[i][0]
        Mi = self.monodromy_matrices[i]
        vj = self.thimbles[j][0]
        Mj = self.monodromy_matrices[j]
        
        res = ((Mi-1)*vi)*self.fiber.intersection_product*(Mj-1)*vj
        
        return 0 if i==j else res*(1 if i>j else -1)
    
    def _restrict_form(self, A):
        """ Given a form A, returns the form A_t such that A/P^k w_n = A_t/P_t^k w_{n-1}dt
        """
        assert self.dim !=0, "cannot restrict form of a dimension 0 variety"
        dt = self._RtoS(self._vars[-1]).derivative() # this formula is specific to Lefschetz fibration where one of the maps is z and the other map does not depend on z, need to change it
        return self._RtoS(A)*dt

    @classmethod
    def _derivative(self, A, P): 
        """computes the numerator of the derivative of A/P^k"""
        field = P.parent()
        return field(A).derivative() - A*P.derivative()         
    
    @classmethod
    def _residue_form(self, A, P, k, alphas): 
        """ returns the formal residue of A/P^k at alpha for alpha in alphas """
        G,U,V = xgcd(P, P.derivative())
        assert G==1
        if k==1:
            return [V(alpha)*A(alpha) for alpha in alphas]
        return _residue_form(A*U/(k-1)+(A*V).derivative()/(k-1)**2, P, k-1, alphas)
    
    @parallel(4)
    @classmethod
    def _compute_transition_matrix_delaunay(cls, L, l, values, nbits=400):
        """computes the numerical transition matrix of L along l, adapted to computations of Delaunay. Accepts l=[]"""
        res = L.numerical_transition_matrix([values[v] for v in l], eps=2**(-nbits), assume_analytic=True) if l!= [] else identity_matrix(L.order())
        return res

    @parallel(4)
    @classmethod
    def _compute_transition_matrix_voronoi(cls, i, L, l, nbits=300):
        """ Returns the numerical transition matrix of L along l, adapted to computations of Voronoi. Accepts l=[]
        """
        logger.info("[%d] Starting integration along edge [%d/%d]"% (os.getpid(), i[0]+1,i[1]))
        res = L.numerical_transition_matrix(l, eps=2**(-nbits), assume_analytic=True) if l!= [] else identity_matrix(L.order())
        logger.info("[%d] Finished integration along edge [%d/%d]"% (os.getpid(), i[0]+1,i[1]))
        return res
               
    def _compute_picard_fuchs(self):
        """ Returns the Picard-Fuchs equations of image of the forms of the variety in the parametrized family of fibers
        """
        # It is arguably more efficient to compute all the Picard-Fuchs equations at the same time
        logger.info("Computing Picard-Fuchs equations of %d forms in dimension %d"% (len(self.cohomology.basis()), self.dim))
        coordinates, denom = self.family.coordinates([self._restrict_form(w) for w in self.cohomology.basis()])
        Ls = [self.family.picard_fuchs_equation(v)*denom for v in coordinates.rows()]
        return Ls
   

    def _compute_paths(self):
        if self.ctx.method=="delaunay":
            self._sps = self._compute_paths_delaunay(self.critical_points)
            paths=[sp+["Rd"]+list(reversed(sp)) for sp in self._sps]
        elif self.ctx.method=="voronoi":
            self._sps, self._loops = self._compute_paths_voronoi([self.ctx.CF(c) for c in self.critical_points]) #basepoint is chosen here
            paths=[list(self._sps[i]+self._loops[i][1:]+list(reversed(self._sps[i]))) for i in range(len(self.critical_points))]
        return paths
 
    def _compute_paths_delaunay(self, singus, shift=1):
        logger.info("Computing paths for integration")
        
        # Need to have a better choice for basepoint, often it is quite far away from other singularities 
        # (on the other hand this sort of guarantees that the degree of basepoint in the cover tree is 1, 
        # which is good for the circuit)
        self.basepoint = floor(min([self.ctx.CF(s).real() for s in singus])-shift) 
        points = [self.basepoint] + singus

        logger.info("Computing Delaunay triangulation of %d points"% (len(singus)+1))
        D = Edges.delaunay(points)
        logger.info("Computing minimal cover tree")
        T = Edges.minimal_cover_tree(D, points)
        logger.info("Computing circuit")
        sps, order = Edges.circuit_delaunay(T, list(range(1, len(points))), 0, points)
        assert len(order)==len(self.critical_points), "circuit returned less critical points than expected"
        self.critical_points = [self.critical_points[i] for i in order]
        orderi = Util.invert_permutation([0] + [i+1 for i in order])
        sps = [[orderi[p] for p in sps[i]] for i in order]
        self.edges = [(orderi[e[0]], orderi[e[1]]) for e in T.edges(labels=False)]
        
        logger.info("Paths for integration are computed")

        return sps
    
    def _compute_paths_voronoi(self,singus, shift=1):

        self.basepoint = floor(min([s.real() for s in singus])-shift)

        logger.info("Computing homotopy representants of the image of the projection in dimension %d"%  self.dim)
        logger.info("Computing Voronoi Diagram of %d points"% len(singus))

        voronoi = Edges.voronoi(singus, self.basepoint)

        sps, loops, order = Edges.voronoi_loops(singus, self.basepoint)
        
        self.critical_points = [self.critical_points[i] for i in order] 
        
        G= Graph(loops=True)
        for sp in sps:
            for i in range(len(sp)-1):
                G.add_edge((sp[i], sp[i+1]))
        for loop in loops:
            for i in range(len(loop)-1):
                G.add_edge((loop[i], loop[i+1]))
            G.add_edge((loop[-1], loop[0]))
        G.remove_loops()
        self.edges = G.edges(labels=False)
        
        return [sps[i] for i in order], [loops[i] for i in order]


