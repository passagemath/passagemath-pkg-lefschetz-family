# `Fibration` documentation

This class allows to compute the monodromy representation of one parameter families of generically smooth projective hypersurfaces.

The first step is to define the polynomial $P$ defining the projective hypersurface $X=V(P)$. For instance, the following gives the Fermat elliptic curve:
```python
R.<X,Y,Z> = PolynomialRing(QQ)
S.<t> = PolynomialRing(R)
P = X**4+Y**4+Z**4 + t*X*Y*Z**2
```
Then the following creates an object representing the fibration:
```python
from lefschetz_family import Fibration
X = Fibration(P)
```
The monodromy representation of the family defined by P is obtained by:
```python
X.monodromy_matrices
```

The module automatically uses available cores for computing numerical integrations and braids of roots. For this, the sage session needs to be made aware of the available cores. This can be done by adding the following line of code before launching the computation (replace `10` by the number of cores you want to use).
```python
os.environ["SAGE_NUM_THREADS"] = '10'
```


## Options
The object `Fibration` can be called with several options:
- `nbits` (positive integer, `400` by default): the number of bits of precision used as input for the computations. If a computation fails to recover the 
- `basepoint` (rational number): which point to use as a basepoint.
- `fibre` (`Hypersurface` object, see above): if the periods of the fiber have already been computed, passing this argument avoids recomputing them. This can only be set if `basepoint` is also set.
- `family` (`Family` object): similarly to `fibre`, one can pass down a `Family` object to avoid redundant computations of Picard-Fuchs equations.
- `cyclic_forms` (list of vectors defined over `self.family.upolring`): a list of cohomology classes of `self.family` that will be used to compute the monodromy. This family has to generate the full cohomology (as a $D$-module).
- `fibration` (list of vectors defined over the rationals): the fibration to pass down to the fibre for computing its periods, see `Hypersurface` above.

## Properties
The object `Fibration` has several properties.
Fibration related properties:
- `critical_values`: the list critical values  of that map.
- `basepoint`: the basepoint of the fibration (i.e. a non critical value).
- `fundamental_group`: the class computing representants of the fundamental group of $\mathbb P^1$ punctured at the critical values.
- `paths`: the list of simple loops around each point of `critical_values`. When this is called, the ordering of `critical_values` changes so that the composition of these loops is the loop around infinity.
- `family`: the one parameter family corresponding to the fibration.
- `fibre`: the fibre above the basepoint.

Miscellaneous properties:
- `P`: the defining equation of $X$.
- `ctx`: the options of $X$, see related section above.

## Copy-paste ready examples

### A family of quartic curves

```python
os.environ["SAGE_NUM_THREADS"] = '10'
from lefschetz_family import Fibration
R.<X,Y,Z> = QQ[]
S.<t> = R[]
P = X**4 + Y**4 + Z**4 + t*X*Y*Z**2
fib = Fibration(P, fibration=[vector([2,0,1]), vector([0,1,0])])
fib.monodromy_matrices
```

### Symmetric cubic surfaces

```python
os.environ["SAGE_NUM_THREADS"] = '10'
from lefschetz_family import Fibration
R.<X,Y,Z,W> = QQ[]
S.<t> = R[]
P = X**3 + Y**3 + Z**3 + W**3 + t*(X*Y*Z + X*Y*W + X*Z*W + Y*Z*W)
fibration = [vector([9, 5, 2, 8]), vector([-7, -9, -8, -1]), vector([6, 6, 9, -8])]
fib = Fibration(P, nbits=800, fibration=fibration)
MatrixGroup(fib.monodromy_matrices).group_id() # we recover the Klein four-group in a few minutes
```
