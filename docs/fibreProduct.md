# `FibreProduct` documentation

This class allows to compute the periods of fibre products of elliptic surfaces.

The first step is to define the elliptic surfaces $S_1$ and $S_2$ defining the fibre product $X=S_1\times_{\mathbb P^1}S_2$ using `EllipticSurface`. It is necessary to give the two elliptic surfaces the same basepoint. The following constructs the example $A\times_{1}c$ of [Periods of fibre products of elliptic surfaces and the Gamma conjecture](https://arxiv.org/abs/2505.07685), Section 6:
```python
from lefschetz_family import EllipticSurface
R.<X,Y,Z> = PolynomialRing(QQ)
S.<t> = R[]
basepoint = -5
S1 = EllipticSurface((-X^2*Z - Y*Z^2)*t - X^3 - X*Y*Z + Y^2*Z, basepoint=basepoint)
S1 = EllipticSurface(Z^3*t^6 - 3*X*Z^2*t^4 + 2*Y*Z^2*t^3 + 3*X^2*Z*t^2 - 3*X*Y*Z*t - X^3 + X*Y*Z + Y^2*Z, basepoint=basepoint)
```
Then the following creates an object representing the hypersurface:
```python
from lefschetz_family import FibreProduct
X = FibreProduct(S1, S2)
```
The period matrix of $X$ is then simply given by:
```python
X.period_matrix
```

The module automatically uses available cores for computing numerical integrations and braids of roots. For this, the sage session needs to be made aware of the available cores. This can be done by adding the following line of code before launching the computation (replace `10` by the number of cores you want to use).
```python
os.environ["SAGE_NUM_THREADS"] = '10'
```


## Copy-paste ready examples

### The Hadamard product $A\times_1 c$
```python
os.environ["SAGE_NUM_THREADS"] = '10'
from lefschetz_family import EllipticSurface
from lefschetz_family import FibreProduct
R.<X,Y,Z> = PolynomialRing(QQ)
S.<t> = R[]
basepoint = -5
S1 = EllipticSurface((-X^2*Z - Y*Z^2)*t - X^3 - X*Y*Z + Y^2*Z, basepoint=basepoint, fibration=[vector(ZZ, [9, -6, -8]), vector(ZZ, [-9, -8, -3])])
S1 = EllipticSurface(Z^3*t^6 - 3*X*Z^2*t^4 + 2*Y*Z^2*t^3 + 3*X^2*Z*t^2 - 3*X*Y*Z*t - X^3 + X*Y*Z + Y^2*Z, basepoint=basepoint, fibration = [vector(ZZ, [-6, -1, -5]), vector(ZZ, [2, -3, 7])])
X = FibreProduct(S1, S2, nbits=1500)
X.period_matrix
```


## Options
The object `FibreProduct` can be called with several options:
- `nbits` (positive integer, `400` by default): the number of bits of precision used as input for the computations. If a computation fails to recover the integral  monodromy matrices, you should try to increase this precision. The output precision seems to be roughly linear with respect to the input precision.
- `debug` (boolean, `False` by default): whether coherence checks should be done earlier rather than late. We recommend setting to true only if the computation failed in normal mode.
- `method` (`"voronoi"` by default/`"delaunay"`/`"delaunay_dual"`): the method used for computing a basis of homotopy. `voronoi` uses integration along paths in the voronoi graph of the critical points; `delaunay` uses integration along paths along the delaunay triangulation of the critical points; `delaunay_dual` paths are along the segments connecting the barycenter of a triangle of the Delaunay triangulation to the middle of one of its edges. In practice, `delaunay` is more efficient for low dimension and low order varieties (such as degree 3 curves and surfaces, and degree 4 curves). This gain in performance is however hindered in higher dimensions because of the algebraic complexity of the critical points (which are defined as roots of high order polynomials, with very large integer coefficients). <b>`"delaunay"` method is not working for now</b>

## Properties
The object `Hypersurface` has several properties.
Fibration related properties, in positive dimension:
- `critical_values`: the list critical values  of that map.
- `basepoint`: the basepoint of the fibration (i.e. a non critical value).
- `fundamental_group`: the class computing representants of the fundamental group of $\mathbb P^1$ punctured at the critical values.
- `paths`: the list of simple loops around each point of `critical_values`. When this is called, the ordering of `critical_values` changes so that the composition of these loops is the loop around infinity.
- `family`: the one parameter family corresponding to the fibration.

Homology related properties:
- `monodromy_matrices`: the matrices of the monodromy action of `paths` on $H_{2}(F_b)$.
- `vanishing_cycles`: the vanshing cycles at each point of `critical_values` along `paths`.
- `thimbles`: the thimbles of $H_3(T^*,F_b)$. They are represented by a starting cycle in $H_2(F_b)$ and (the index of) an element of `paths`.
- `extensions`: integer linear combinations of thimbles with vanishing boundary.
- `infinity_loops`: extensions around the loop at infinity.
- `homology`: a basis of $\Lambda_{\rm vc}^\perp/\Lambda_{\rm vc}$.
- `intersection_product`: the intersection product of $H_3(\Lambda_{\rm vc}^\perp/\Lambda_{\rm vc})$.
- `lift`: a map taking a linear combination of thimbles with zero boundary (i.e. an element of $\ker\left(\delta:H_3(T^*, F_b)\to H_{2}(F_b)\right)$) and returning the homology class of its lift in $H_3(\Lambda_{\rm vc}^\perp/\Lambda_{\rm vc})$, in the basis `homology_modification`.
- `types`: `types[i]` is the type of the fibre above `critical_values[i]`.

Cohomology related properties:
- `cohomology`: the elements of the cohomology in which `period_matrix` is given, represented by triples `w1, w2, f` representing the form $\frac{1}{f}\omega_1\otimes\omega_2\wedge {\mathrm dt}$. There is no garantee that this yields a basis.
- `picard_fuchs_equations`: the Picard-Fuchs equations of the elements of `cohomology`

Period related properties
- `period_matrix`: the period matrix of $X$ in the aforementioned bases `homology` and `cohomology`, as well as the cohomology class of the linear section in even dimension

Miscellaneous properties:
- `S1` and `S2`: the underlying elliptic surfaces.
- `dim`: the dimension of $X$.
- `degree`: the degree of $X$.
- `ctx`: the options of $X$, see related section above.
