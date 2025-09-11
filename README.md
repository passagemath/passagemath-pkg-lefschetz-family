# lefschetz-family


## Description
This Sage package provides a means of efficiently computing periods of complex projective hypersurfaces and elliptic surfaces over $\mathbb P^1$ with certified rigorous precision bounds.
It implements the methods described in 
- [Effective homology and periods of complex projective hypersurfaces](https://doi.org/10.1090/mcom/3947) ([arxiv:2306.05263](https://doi.org/10.48550/arXiv.2306.05263)).
- [A semi-numerical algorithm for the homology lattice and periods of complex elliptic surfaces over the projective line](https://doi.org/10.1016/j.jsc.2024.102357) ([arxiv:2401.05131](https://arxiv.org/abs/2401.05131)).
- [Periods of fibre products of elliptic surfaces and the Gamma conjecture](https://doi.org/10.48550/arXiv.2505.07685) ([arxiv:2505.07685](https://doi.org/10.48550/arXiv.2505.07685)).
- [Periods in algebraic geometry : computations and application to Feynman integrals](https://theses.hal.science/tel-04823423) ([hal:tel-04823423](https://theses.hal.science/tel-04823423)).
- [Galois groups of symmetric cubic surfaces](https://doi.org/10.48550/arXiv.2509.06785) ([arxiv:2509.06785](https://doi.org/10.48550/arXiv.2509.06785)).

Please cite accordingly.

This package is a successor to the [numperiods](https://gitlab.inria.fr/lairez/numperiods) package by Pierre Lairez. It contains files taken from this package, that have sometimes been slightly modified to accomodate for new usage.

## How to install

In a terminal, run
```
sage -pip install git+https://github.com/mkauers/ore_algebra.git
sage -pip install lefschetz-family
```
or
```
sage -pip install --user git+https://github.com/mkauers/ore_algebra.git
sage -pip install --user lefschetz-family
```

Alternatively, install the `ore_alegbra` package (available at [https://github.com/mkauers/ore_algebra](https://github.com/mkauers/ore_algebra)), then download this repository and add the path to the main folder to your `sys.path`.

## Requirements
Sage 9.0 and above is recommended. Furthermore, this package has the following dependencies:

- [Ore Algebra](https://github.com/mkauers/ore_algebra).
- The [delaunay-triangulation](https://pypi.org/project/delaunay-triangulation/) package from PyPI.



## Documentation

- [Hypersurface](docs/hypersurface.md) for computing periods of hypersurfaces.
- [EllipticSurface](docs/ellipticSurface.md) for computing  period of elliptic surfaces.
- [DoubleCover](docs/doubleCover.md) for computing periods of ramified double cover of projective spaces.
- [FibreProduct](docs/fibreProduct.md) for computing periods of fibre products of elliptic surfaces.
- [Fibration](docs/fibration.md) for computing monodromy representations of families of hypersurfaces.

## Performance benchmarking

Here is a runtime benchmark for various examples, with an input precision of 1000 bits:
| Hypersurface (generic) 	| Time (on 10 M1 cores) | Recovered precision (decimal digits)  |
|-------------------	|----------------------	| ----------------------                |
| Elliptic curve    	| 5 seconds            | 340 digits                            |
| Quartic curve     	| 90 seconds            	| 340 digits                            |
| Quintic curve     	| 5 minutes            	| 340 digits                            |
| Sextic curve      	| 30 minutes            	| 300 digits                            |
| Cubic surface     	| 40 seconds         	| 340 digits                            |
| Quartic surface   	| 1 hour        	    | 300 digits                            |
| Cubic threefold   	| 15 minutes        	    | 300 digits                            |
| Cubic fourfold   	| 20 hours        	    | 300 digits                            |
| Rational elliptic surface | 10 seconds        	    | 300 digits                            |
| Elliptic K3 surface   	| 30 seconds*        	    | 300 digits                            |
| Degree 2 K3 surface   	| 5 minutes        	    | 300 digits                            |


*for holomorphic periods



## Contact
For any question, bug or remark, please contact [eric.pichon@mis.mpg.de](mailto:eric.pichon@mis.mpg.de).

## Roadmap
Near future milestones:
- [x] Encapsulate integration step in its own class
- [x] Certified computation of the exceptional divisors
- [x] Saving time on differential operator by precomputing cache before parallelization
- [x] Computing periods of elliptic fibrations.
- [x] Removing dependency on `numperiods`.

Middle term goals include:
- [ ] Making Delaunay triangulation functional again
- [ ] Having own implementation of 2D voronoi graphs/Delaunay triangulation

Long term goals include:
- [x] Tackling cubic threefolds.
- [x] Generic code for all dimensions.
- [x] Computing periods of K3 surfaces with mildy singular quartic models.
- [ ] Dealing with other singularities, especially curves.
- [ ] Computing periods of complete intersections.
- [x] Computing periods of weighted projective hypersurfaces, notably double covers of $\mathbb P^2$ ramified along a sextic.

Other directions include:
- [ ] Computation of homology through braid groups instead of monodromy of differential operators.


## Project status
This project is actively being developped.
