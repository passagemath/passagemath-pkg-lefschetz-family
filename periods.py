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
from voronoi import FundamentalGroupVoronoi
from integrator import Integrator
from Util import Util
from Context import Context

import logging
import os
import time

logger = logging.getLogger(__name__)


class LefschetzFamily(object):
    def __init__(self, P, axis=None, **kwds):
        """P, a homogeneous polynomial defining a smooth hypersurface X in P^{n+1}.

        This class aims at computing an effective basis of the homology group H_n(X), 
        given as lifts of paths through a Lefschetz fibration.
        """
        
        self.ctx = Context(**kwds)
        
        # assert P.is_homogeneous(), "nonhomogeneous defining polynomial"
        
        self._P = P
        self._axis=axis
    
    
    @property
    def intersection_product(self):
        if not hasattr(self,'_intersection_product'):
            # assert self.dim==0, "No way to compute intersection product in positive dimension yet"
            if self.dim==0:
                self._intersection_product=identity_matrix(self.degree)
            else:
                self._intersection_product=self._compute_intersection_product()
        return self._intersection_product

    @property
    def period_matrix(self):
        if not hasattr(self, '_period_matrix'):
            if self.dim==0:
                R = self.P.parent()
                affineR = PolynomialRing(QQbar, 'X')
                affineProjection = R.hom([affineR.gens()[0],1], affineR)
                self._period_matrix = matrix([self._residue_form(affineProjection(b), affineProjection(self.P), (b.degree()+len(R.gens()))//self.P.degree(), self.homology) for b in self.cohomology]).change_ring(self.ctx.CBF)
            else:
                self._period_matrix = matrix(self.integrated_thimbles([i for i in range(len(self.cohomology))]))*matrix(self.homology).transpose()

        return self._period_matrix

    @property
    def simple_periods(self):
        if not hasattr(self, '_simple_periods'):
            if self.dim==0:
                self._simple_periods = self.period_matrix
            else:
                self._simple_periods = matrix(self.integrated_thimbles([0]))*matrix(self.homology).transpose()
        return self._simple_periods



    @property
    def P(self):
        return self._P

    @property
    def degree(self):
        if not hasattr(self,'_degree'):
            self._degree = self.P.degree()
        return self._degree
    
    @property
    def dim(self):
        if not hasattr(self,'_dim'):
            self._dim = len(self.P.parent().gens())-2
        return self._dim

    def picard_fuchs_equation(self, i):
        if not hasattr(self,'_picard_fuchs_equations'):
            _picard_fuchs_equations = [None for i in range(len(self.cohomology))]
            logger.info("Computing Picard-Fuchs equations of %d forms in dimension %d"% (len(self.cohomology), self.dim))
            coordinates, denom = self.family.coordinates([self._restrict_form(w) for w in self.cohomology])

            for j, v in enumerate(coordinates.rows()):
                v2 = v/denom
                denom2 = lcm([r.denominator() for r in v2 if r!=0])
                numerators = denom2 * v2
                L = self.family.picard_fuchs_equation(numerators)*denom2
                L = DifferentialOperator(L)
                logger.info("Operator [%d/%d] has order %d and degree %d for form with numerator of degree %d"% (j+1, len(self.cohomology), L.order(), L.degree(), self.cohomology[j].degree()))
                _picard_fuchs_equations[j] = L
                self._picard_fuchs_equations = _picard_fuchs_equations
        return self._picard_fuchs_equations[i]
    
    @property
    def cohomology(self):
        if not hasattr(self,'_cohomology'):
            self._cohomology = Cohomology(self.P).basis()
        return self._cohomology
    
    
    @property
    def family(self):
        if not hasattr(self,'_family'):
            RtoS = self._RtoS()
            self._family = Family(RtoS(self.P))
        return self._family
    

    @property
    def fibration(self):
        if not hasattr(self,'_fibration'): #TODO try to reduce variance of distance between critical points(?)
            l = vector([randint(-10,10) for i in range(self.dim+2)]) if self._axis==None else self._axis
            m = vector([randint(-10,10) for i in range(self.dim+2)])
            assert matrix([l,m]).rank()==2, "fibration is not well defined"
            self._fibration= (l,m)
        return self._fibration

    # def _compute_distance(self, points):
    #     return min([min([abs(ComplexField(50)(p1-p2)) for p2 in points if p2!=p1]) for p1 in points])
    

    @property
    def critical_points(self):
        if not hasattr(self,'_critical_points'):
            R = self.P.parent()
            _vars = [v for v in R.gens()]
            logger.info("Computing critical points")
            begin=time.time()
            forms=[v.dot_product(vector(_vars)) for v in self.fibration]
            f=forms[0]/forms[1]
            S = PolynomialRing(QQ, _vars+['k','t'])
            k,t= S.gens()[-2:]
            eqs = [
                self.P, 
                forms[1]-1, 
                t*forms[1]-forms[0]
            ] + [(f.derivative(var).numerator()*k-self.P.derivative(var)*f.derivative(var).denominator()) for var in _vars]

            ideal = S.ideal(eqs).elimination_ideal(S.gens()[:-1])
            Qt = PolynomialRing(QQ, 't')

            roots_with_multiplicity = Qt(ideal.groebner_basis()[0]).roots(AlgebraicField())
            if not self.ctx.debug:
                for e in roots_with_multiplicity:
                    assert e[1]==1, "double critical values, fibration is not Lefschetz"
            self._critical_points=[e[0] for e in roots_with_multiplicity]
            end = time.time()
            logger.info("Critical points computed in %s"% time.strftime("%H:%M:%S",time.gmtime(end-begin)))
            # if self.dim==2:
            #     assert len(self._critical_points)==36, "fibration is not Lefschetz (if you are not dealing with quartics, this should be removed)"
        return self._critical_points
    
    @property
    def monodromy_matrices(self):
        assert self.dim!=0, "Cannot compute monodromy matrices in dimension 0"
        if not hasattr(self, '_monodromy_matrices'):
            i=0
            assert self.picard_fuchs_equation(i).order()== len(self.family.basis),"Picard-Fuchs equation is not cyclic, cannot use it to compute monodromy"
            transition_matrices= self.transition_matrices([i])[0]

            n = len(self.fiber.homology)
            r = len(self.critical_points)  

            if self.dim%2==1:
                invariant_vector=vector([1]*n)
                proj = block_matrix([[identity_matrix(n-1), matrix([[-1]]*(n-1))]],subdivide=False)
                lift = block_matrix([[identity_matrix(n-1)], [matrix([[0]*(n-1)])]],subdivide=False)
                basis_change = block_matrix([[lift, matrix([-invariant_vector]).transpose()]],subdivide=False)
            
            logger.info("Computing the coordinates of the successive derivatives of integration forms")
            begin = time.time()
            # derivatives_coordinates, denom = self.derivatives_coordinates(i)
            derivatives_at_basepoint = self.derivatives_values_at_basepoint(i)
            end = time.time()
            duration_str = time.strftime("%H:%M:%S",time.gmtime(end-begin))
            logger.info("Coordinates computed in %s"% duration_str)
            
            integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(n+1 if self.dim%2==0 else n)])
            # initial_conditions = integration_correction* derivatives_coordinates(self.basepoint)/denom(self.basepoint)*self.fiber.period_matrix
            initial_conditions = integration_correction* derivatives_at_basepoint*self.fiber.period_matrix

            if self.dim%2==1:
                initial_conditions = initial_conditions*lift

            Ms = [(initial_conditions.submatrix(1,0)**(-1)*M.submatrix(1,1)*initial_conditions.submatrix(1,0)) for M in transition_matrices]
            if not self.ctx.debug:
                Ms = [M.change_ring(ZZ) for M in Ms]

            for i in range(r):
                M=Ms[i]
                if self.dim%2==1:
                    prod= invariant_vector*self.fiber.intersection_product*(lift*M*proj-identity_matrix(n))*lift/n
                    Ms[i] = basis_change.inverse()*block_matrix([[M,zero_matrix(n-1,1)],[matrix([prod]),matrix([[1]])]])*basis_change

                if not self.ctx.singular and not self.ctx.debug:
                    assert (Ms[i]-1).rank()==1, "If M is a monodromy matrix around a single critical point, M-I should have rank 1"
            self._monodromy_matrices = Ms
        return self._monodromy_matrices
    
    @property
    def fiber(self):
        assert self.dim!=0, "Variety of dimension 0 does not have a fiber"
        if not hasattr(self,'_fiber'):
            RtoS = self._RtoS()
            evaluate_at_basepoint = RtoS.codomain().hom([self.basepoint], RtoS.codomain().base_ring())
            self._fiber = LefschetzFamily(evaluate_at_basepoint(RtoS(self.P)), method=self.ctx.method, nbits=self.ctx.nbits, depth=self.ctx.depth+1)

        return self._fiber

    @property
    def thimbles(self):
        if not hasattr(self,'_thimbles'):
            self._thimbles=[]
            for pc, path in zip(self.permuting_cycles, self.paths):
                self._thimbles+=[(pc, path)]
        return self._thimbles

    @property
    def permuting_cycles(self):
        if not hasattr(self, '_permuting_cycles'):
            self._permuting_cycles = [None for i in range(len(self.monodromy_matrices))]
            for i in range(len(self.monodromy_matrices)):
                M = self.monodromy_matrices[i]
                D, U, V = (M-1).smith_form()
                self._permuting_cycles[i] = V * vector([1]+[0]*(V.dimensions()[0]-1))
        return self._permuting_cycles

    @property
    def vanishing_cycles(self):
        if not hasattr(self, '_vanishing_cycles'):
            self._vanishing_cycles = []
            for p, M in zip(self.permuting_cycles,self.monodromy_matrices):
                self._vanishing_cycles += [(M-1)*p]
        return self._vanishing_cycles


    @property
    def infinity_loops(self):
        if not hasattr(self, '_infinity_loops'):
            Mtot=1
            phi=[]
            for M, v in zip(self.monodromy_matrices, self.vanishing_cycles):
                tempM=(M-1)*Mtot
                phi+=[[c/v for c in tempM.columns()]]
                Mtot=M*Mtot
            phi = matrix(phi).transpose().change_ring(ZZ)
            if not self.ctx.debug:
                assert Mtot == identity_matrix(len(self.fiber.homology)), "Monodromy around infinity is nontrivial, most likely because the paths do not actually compose to the loop around infinity"
            self._infinity_loops = phi.rows()

        return self._infinity_loops
    
    @property
    def extensions(self):
        if not hasattr(self, '_extensions'):
            delta = matrix(self.vanishing_cycles).change_ring(ZZ)
            self._extensions = delta.kernel()

        return self._extensions
    

    @property
    def homology(self):
        if not hasattr(self, '_homology'):
            if self.dim==0:
                R = self.P.parent()
                affineR = PolynomialRing(QQbar, 'X')
                affineProjection = R.hom([affineR.gens()[0],1], affineR)
                self._homology = [e[0] for e in affineProjection(self.P).roots()]

            else:
                r = len(self.monodromy_matrices)
                n = len(self.fiber.cohomology)

                # compute representants of the quotient H(Y)/kerdelta
                D, U, V = self.extensions.matrix().smith_form()
                B = D.solve_left(matrix(self.infinity_loops)*V).change_ring(ZZ)*U
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
                self._homology = (quotient_basis*self.extensions.matrix()).rows() # NB this is the homology of Y, to recover the homology of X we need to remove the kernel of the period matrix
        return self._homology

    @property
    def exceptional_divisors(self):
        if not hasattr(self, '_exceptional_divisors'):
            assert self.dim<=2, "Not implemented yet"
            self._exceptional_divisors = []
        return self._exceptional_divisors


    def _RtoS(self):
        R = self.P.parent()
        for varid in range(self.dim+2):
            if self.fibration[0][varid] !=0:
                break
        self._restriction_variable = varid # this is the variable we're removing when going from P to Pt

        a = self.fibration[1][self._restriction_variable]
        b = self.fibration[0][self._restriction_variable]

        _vars = [v for v in R.gens()]
        S = PolynomialRing(PolynomialRing(QQ, [_vars[i] for i in range(len(_vars)) if i != self._restriction_variable]), 't')
        t=S.gens()[0]

        l = vector([_vars[i] for i in range(len(_vars)) if i != self._restriction_variable])*vector([self.fibration[0][i] for i in range(len(_vars)) if i != self._restriction_variable])
        m = vector([_vars[i] for i in range(len(_vars)) if i != self._restriction_variable])*vector([self.fibration[1][i] for i in range(len(_vars)) if i != self._restriction_variable])

        form = (-S(l)+t*S(m))
        denom = b-a*t

        RtoS = R.hom([denom*S(_vars[i]) if i != self._restriction_variable else form for i in range(len(_vars))], S)

        return RtoS

    def _restrict_form(self, A):
        """ Given a form A, returns the form A_t such that A/P^k w_n = A_t/P_t^k w_{n-1}dt
        """
        assert self.dim !=0, "cannot restrict form of a dimension 0 variety"

        RtoS = self._RtoS()


        R = self.P.parent()

        a = self.fibration[1][self._restriction_variable]
        b = self.fibration[0][self._restriction_variable]

        _vars = [v for v in R.gens()]
        S = PolynomialRing(PolynomialRing(QQ, [_vars[i] for i in range(len(_vars)) if i != self._restriction_variable]), 't')
        t=S.gens()[0]

        l = vector([_vars[i] for i in range(len(_vars)) if i != self._restriction_variable])*vector([self.fibration[0][i] for i in range(len(_vars)) if i != self._restriction_variable])
        m = vector([_vars[i] for i in range(len(_vars)) if i != self._restriction_variable])*vector([self.fibration[1][i] for i in range(len(_vars)) if i != self._restriction_variable])

        form = (-S(l)+t*S(m))
        denom = b-a*t

        dt = (-a * S(l) +b* S(m))*denom**(self.dim)
        return RtoS(A)*dt

    def transition_matrices(self, l):
        if not hasattr(self, '_integratedQ'):
            self._integratedQ = [False for i in range(len(self.cohomology))]
        if not hasattr(self, '_transition_matrices'):
            self._transition_matrices = [None for i in range(len(self.cohomology))]
        for i in l:
            if not self._integratedQ[i]:
                L = self.picard_fuchs_equation(i)
                L = L* L.parent().gens()[0]
                self._transition_matrices[i] = self.integrate(L)
                self._integratedQ[i]=True
        return [self._transition_matrices[i] for i in l]
    
    def integrate(self, L):
        logger.info("Computing numerical transition matrices of operator of order %d and degree %d (%d edges total)."% (L.order(), L.degree(), len(self.fundamental_group.edges)))
        begin = time.time()

        integrator = Integrator(self.fundamental_group, L, self.ctx.nbits)
        transition_matrices = integrator.transition_matrices
        
        end = time.time()
        duration_str = time.strftime("%H:%M:%S",time.gmtime(end-begin))
        logger.info("Integration finished -- total time: %s."% duration_str)

        return transition_matrices

    def forget_transition_matrices(self):
        del self._integratedQ
        del self._transition_matrices
        del self._integrated_thimbles
        del self._integrated_thimblesQ

    def integrated_thimbles(self, l):
        transition_matrices= self.transition_matrices(l)
        if not hasattr(self, '_integrated_thimblesQ'):
            self._integrated_thimblesQ = [False for i in range(len(self.cohomology))]
        if not hasattr(self, '_integrated_thimbles'):
            self._integrated_thimbles = [None for i in range(len(self.cohomology))]
        
        N=len(self.cohomology)
        n=len(self.fiber.homology)
        r=len(self.thimbles)

        for i2 in range(len(l)):
            i= l[i2]
            if not self._integrated_thimblesQ[i]:
                derivatives_at_basepoint = self.derivatives_values_at_basepoint(i)
                integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(n+1 if self.dim%2==0 else n)])
                initial_conditions = integration_correction* derivatives_at_basepoint*self.fiber.period_matrix
                self._integrated_thimbles[i]=[(transition_matrices[i2][j]*initial_conditions*self.permuting_cycles[j])[0] for j in range(r)]
                self._integrated_thimblesQ[i] = True
        return [self._integrated_thimbles[i] for i in l]
    
    # Integration methods

    def derivatives_coordinates(self, i):
        if not hasattr(self, '_coordinatesQ'):
            self._coordinatesQ = [False for i in range(len(self.cohomology))]
        if not hasattr(self, '_coordinates'):
            self._coordinates = [None for i in range(len(self.cohomology))]

        if not self._coordinatesQ[i]:
            n=len(self.fiber.homology)
            RtoS = self._RtoS()
            w = self.cohomology[i]
            wt = self._restrict_form(w)
            derivatives = [RtoS(0), wt]
            for k in range(n-1 if self.dim%2==0 else n-2):
                derivatives += [self._derivative(derivatives[-1], RtoS(self.P))] 
            self._coordinates[i] = self.family.coordinates(derivatives)
            self._coordinatesQ[i] = True

        return self._coordinates[i]


    def derivatives_values_at_basepoint(self, i):
        RtoS = self._RtoS()
        n=len(self.fiber.homology)


        w = self.cohomology[i]
        wt = self._restrict_form(w)
        derivatives = [RtoS(0), wt]
        for k in range(n-1 if self.dim%2==0 else n-2):
            derivatives += [self._derivative(derivatives[-1], RtoS(self.P))] 
        return self.family._coordinates(derivatives, self.basepoint)

    def _compute_intersection_product(self):
        r=len(self.thimbles)
        inter_prod_thimbles = matrix([[self._compute_intersection_product_thimbles(i,j) for j in range(r)] for i in range(r)])
        intersection_11 = (-1)**self.dim * (matrix(self.homology)*inter_prod_thimbles*matrix(self.homology).transpose()).change_ring(ZZ)
        if self.dim%2==0:
            intersection_02 = zero_matrix(2,2)
            intersection_02[0,1], intersection_02[1,0] = 1,1
            intersection_02[1,1] = -1
            return block_diagonal_matrix(intersection_11, intersection_02)
        else:
            return intersection_11
        
    def _compute_intersection_product_thimbles(self,i,j):
        vi = self.thimbles[i][0]
        Mi = self.monodromy_matrices[i]
        vj = self.thimbles[j][0]
        Mj = self.monodromy_matrices[j]

        di, dj = ((Mi-1)*vi), (Mj-1)*vj

        
        res = di*self.fiber.intersection_product*dj
        resid = -vi*self.fiber.intersection_product*di

        if i==j:
            return resid
        if i<j:
            return res
        else:
            return 0

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

    @property
    def fundamental_group(self):
        if not hasattr(self,'_fundamental_group'):
            begin = time.time()

            fundamental_group = FundamentalGroupVoronoi(self.critical_points, self.basepoint) # access future delaunay implem here
            fundamental_group.sort_loops()

            end = time.time()
            duration_str = time.strftime("%H:%M:%S",time.gmtime(end-begin))
            logger.info("Fundamental group computed in %s."% duration_str)

            self._critical_points = fundamental_group.points[1:]
            self._fundamental_group = fundamental_group
        return self._fundamental_group

    @property
    def paths(self):
        if not hasattr(self,'_paths'):
            paths = []
            for path in self.fundamental_group.pointed_loops:
                paths += [[self.fundamental_group.vertices[v] for v in path]]
            self._paths= paths
        return self._paths

    @property
    def basepoint(self):
        if  not hasattr(self, '_basepoint'):
            shift = 1
            reals = [self.ctx.CF(c).real() for c in self.critical_points]
            xmin, xmax = min(reals), max(reals)
            self._basepoint = Util.simple_rational(xmin - (xmax-xmin)*shift, (xmax-xmin)/10)
        return self._basepoint

