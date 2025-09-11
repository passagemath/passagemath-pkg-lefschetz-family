# `DoubleCover` documentation

This class allows to compute the periods of double covers of projective space.

## Usage

The defining equation for the double cover should be given as a homogeneous polynomial of even degree. Such a polynomial $P$ represents the double cover $X = V(w^2-P)$.
```python
R.<X,Y,Z> = PolynomialRing(QQ)
P = X^6+Y^6+Z^6
```
Then the following creates an object representing the variety:
```python
from lefschetz_family import DoubleCover
X = DoubleCover(P)
```
## Copy-paste ready examples

TODO

## Options

The options are the same as those for `Hypersurface` (see above).

## Properties

The object `DoubleCover` has several properties.
Fibration related properties, in positive dimension:
<!-- - `fibration`: the two linear maps defining the map $X\dashrightarrow \mathbb P^1$. -->
- `critical_values`: the list critical values  of that map.
- `basepoint`: the basepoint of the fibration (i.e. a non critical value).
- `fibre`: the fibre above the basepoint as a `DoubleCover` object.
- `paths`: the list of simple loops around each point of `critical_points`. When this is called, the ordering of `critical_points` changes so that the composition of these loops is the loop around infinity.
- `family`: the one parameter family corresponding to the fibration.

Homology related properties:
- `extensions`: the extensions of the fibration.
- `homology`: the homology of $X$.
- `fibre_class`: the class of the fibre in `homology`.
- `section`: the class of the zero section in `homology`.
- `intersection_product`: the intersection matrix of the surface in the basis `homology`.
- `lift`: a map taking a combination of thimbles of the morsification with empty boundary and returning its class in `homology`.

Cohomology related properties:
- `holomorphic_forms`: a basis of rational functions $f(t)$ such that $f(t) {Res}\frac{\Omega_2}{P_t}\wedge\mathrm dt$ is a holomorphic form of $S$.

Period related properties:
- `period_matrix`: the holomorphic periods of $X$ in the bases `self.homology` and `self.holomorphic_forms`.
- `effective_periods`: the holomorphic periods $X$ in the bases `self.effective_lattice` and `self.holomorphic_forms`

Sublattices of homology. Unless stated otherwise, lattices are given by the coordinates of a basis of the lattice in the basis `homology`:
- `primary_lattice`: The lattice of effective cycles of $X$, consisting of the concatenation of `extensions`, `singular_components`, `fibre_class` and `section`.

Miscellaneous properties:
- `P`: the defining equation of $X$.
- `dim`: the dimension of $X$.
- `degree`: the degree of $P$.
- `ctx`: the options of $X$, see related section above.
