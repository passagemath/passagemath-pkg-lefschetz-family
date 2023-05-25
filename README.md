# lefschetz-family


## Description
This Sage package provides a means of efficiently computing periods of complex projective hypersurfaces with certified rigorous precision bounds.
It implements the methods described in [Coming soon]().
Here is a runtime benchmark for various examples:
| Variety (generic) 	| Time (on 4 M1 cores) 	|
|-------------------	|----------------------	|
| Elliptic curve    	| 7 seconds            	|
| Quartic curve     	| 4 minutes            	|
| Cubic surface     	| 4 minutes 20         	|
| Quartic surface   	| est. 20 hours        	|

## Requirements
Sage 9.0 and above is recommended. Furthermore, this package has the following dependencies:

- [Ore Algebra](https://github.com/mkauers/ore_algebra).
- [numperiods](https://gitlab.inria.fr/lairez/numperiods).

## Usage
The first step is to define the polynomial $P$ defining the projective hypersurface $X=V(P)$. For instance, the following gives the Fermat elliptic curve:
```python
R.<X,Y,Z> = PolynomialRing(QQ)
P = X**3+Y**3+Z**3
```
Then the following creates an object of the period:
```python
from period import LefschetzFamily
X = LefschetzFamily(P)
```
The period matrix of $X$ is the simply given by:
```python
X.period_matrix
```

See [Coming soon]() for a usage example.


### Options
The object `LefschetzFamily` can be called with several options:
- `method` (`"voronoi"` by default/`"delaunay"`): the method used for computing a basis of homotopy. `voronoi` uses integration along paths in the voronoi graph of the critical points, whereas `delaunay` uses integration along paths along the delaunay triangulation of the critical points. In practice, `delaunay` is more efficient for low dimension and low order varieties (such as degree 3 curves and surfaces, and degree 4 curves). This gain in performance is however hindered in higher dimensions because of the algebraic complexity of the critical points (which are defined as roots of high order polynomials, with very large integer coefficients). <b>`"delaunay"` method is not working for now</b>
- `nbits` (positive integer, `400` by default): the number of bits of precision used as input for the computations. If a computation fails to recover the integral  monodromy matrices, you should try to increase this precision. The output precision seems to be roughly linear with respect to the input precision.
- `debug` (boolean, `False` by default): whether coherence checks should be done earlier rather than late. Set to true only if the computation fails.
- `singular` (boolean, `False` by default): whether the variety is singular. <b>Not implemented yet</b>

### Properties
The object `LefschetzFamily` has several properties.
Fibration related properties, in positive dimension:
- `fibration`: the two linear maps defining the map $X\dashrightarrow \mathbb P^1$.
- `critical_points`: the list critical values  of that map.
- `basepoint`: the basepoint of the fibration (i.e. a non critical value).
- `fiber`: the fiber above the basepoint as a `LefschetzFamily` object.
- `paths`: the list of simple loops around each point of `critical_points`. When this is called, the ordering of `critical_points` changes so that the composition of these loops is the loop around infinity.
- `family`: the one parameter family corresponding to the fibration.

Homology related properties:
- `monodromy_matrices`: the matrices of the monodromy action of `paths` on $H_{n-1}(X_b)$.
- `vanishing_cycles`: the vanshing cycles at each point of `critical_points` along `paths`.
- `thimbles`: the thimbles of $H_n(Y,Y_b)$. They are represented by a starting cycle in $H_n(Y_b)$ and a loop in $\mathbb C$ avoiding `critical_points` and pointed at `basepoint`.
- `extensions`: integer linear combinations of thimbles with vanishing boundary.
- `infinity_loops`: extensions around the loop at infinity.
- `exceptional_divisors`: the exceptional cycles coming from the modification $Y\to X$. <b>(WIP)</b>
- `homology`: a basis of representants of $H_n(X)$, apart from the $H_n(X_b)$ term.

Cohomology related properties:
- `cohomology`: a basis of $PH^n(X)$, represented by the numerators of the rational fractions.
- `picard_fuchs_equation(i)`: the picard fuchs equation of the parametrization of i-th element of `cohomology` by the fibration

Period related properties
- `period_matrix`: the period matrix of the blowup of $X$ in the aforementioned bases `homology` and `cohomology`
- `simple_periods`: the periods of the first element of `cohomology` in the basis `homology`.

Miscellaneous properties:
- `dim`: the dimension of $X$.
- `ctx`: the options of $X$, see related section above.


## Contact
For any question, bug or remark, please contact [eric.pichon@polytechnique.edu](mailto:eric.pichon@polytechnique.edu).

## Roadmap
Near future milestones:
- [x] Encapsulate integration step in its own class
- [ ] Certified computation of the exceptional divisors
- [ ] Making Delaunay triangulation functional again
- [x] Saving time on differential operator by precomputing cache before parallelization

Middle term goals include:
- [ ] Having own implementation of 2D voronoi graphs/Delaunay triangulation

Long term goals include:
- [ ] Computing periods of elliptic fibrations.
- [ ] Tackling higher dimensional varieties (most notably cubics in $\mathbb P^5$).
- [ ] Computing periods of singular varieties.
- [ ] Computing periods of complete intersections.
- [ ] Computing periods of weighted projective hypersurfaces, notably double covers of $\mathbb P^2$ ramified along a cubic.

Other directions include:
- [ ] Computation of homology through braid groups instead of monodromy of differential operators.


## Project status
This project is actively being developped.
