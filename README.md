# lefschetz-family


## Description
This Sage package provides a means of efficiently computing periods of complex projective hypersurfaces with certified rigorous precision bounds. Here is a runtime benchmark for various examples:
| Variety (generic) 	| Time (on 4 M1 cores) 	|
|-------------------	|----------------------	|
| Elliptic curve    	| 7 seconds            	|
| Quartic curve     	| 4 minutes            	|
| Cubic surface     	| 4 minutes 20         	|
| Quartic surface   	| est. 20 hours        	|

## Requirements
Sagee 9.0 is recommended. Furthermore, this project relies on the following packages:

- [Ore Algebra](https://github.com/mkauers/ore_algebra). For better performance, I recommend using the [`fastalgexp` branch in Marc Mezzarobba's fork](https://github.com/mezzarobba/ore_algebra/tree/fastalgexp).
- [The branch `coordinates` in my fork of numperiods](https://gitlab.inria.fr/epichonp/numperiods/-/tree/coordinates).
- The [delaunay-triangulation](https://pypi.org/project/delaunay-triangulation/) package from PyPI.

## Usage
The first step is to defined the polynomial defining the projective hypersurface. For instance, the following gives the Fermat elliptic curve:
```python
R.<X,Y,Z> = PolynomialRing(QQ)
P = X**3+Y**3+Z**3
```
Then the following creates an object of the period:
```python
from period import LefschetzFamily
X = LefschetzFamily(P)
```
The period matrix of X is the simply given by:
```python
X.period_matrix
```

### Options
The object `LefschetzFamily` can be called with several options:
- `compute_periods` (`True` by default): whether the code should compute the periods of the variety. If `False`, the algorithm stops once it has computed a basis of the homology of the variety.
- `singular` (`False` by default): Not implemented yet -- whether the variety is singular. If `True`, the computation of cohomology is a bit more costly, and singular critical points do not behave like regular critical points.
- `method` (`voronoi` by default/`delaunay`): the method used for computing a basis of homotopy. `voronoi` uses integration along paths in the voronoi graph of the critical points, whereas `delaunay` uses integration along paths along the delaunay triangulation of the critical points. In practice, `delaunay` is more efficient for low dimension and low order varieties (such as degree 3 curves and surfaces, and degree 4 curves). This gain in performance is however hindered in higher dimensions because of the algebraic complexity of the critical points (which are defined as roots of high order polynomials, with very large integer coefficients).
- `nbits` (`400` by default): the number of bits of precision used as input for the computations. If a computation fails to recover the integral  monodromy matrices, you should try to increase this precision. The output precision seems to be roughly linear with respect to the input precision.

## Contact
For any question, bug or remark, please contact [eric.pichon@polytechnique.edu](mailto:eric.pichon@polytechnique.edu).

## Roadmap
Near future milestones:
- [ ] Encapsulation integration step in its own class
- [ ] Making Delaunay triangulation functional again (issues with graphs)
- [ ] Saving time on differential operator by precomputing cache before parallelization

Middle term goals include:
- [ ] Having own implementation of 2D voronoi graphs/Delaunay triangulation

Long term goals include:
- [ ] Tackling higher dimensional varieties (most notably cubics in P^5).
- [ ] Computing periods of singular varieties.
- [ ] Computing periods of elliptic fibrations.
- [ ] Computing periods of complete intersections.
- [ ] Computing periods of weighted projective hypersurfaces, notably double covers of P^2 ramified along a cubic.

Other directions include:
- [ ] Computation of homology through braid groups instead of monodromy of differential operators


## Project status
This project is actively being developped.
