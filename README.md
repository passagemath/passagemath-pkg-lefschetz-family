# lefschetz-family



## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.inria.fr/epichonp/lefschetz-family/-/settings/integrations)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Description
This Sage package provides a means of efficiently computing periods of complex projective hypersurfaces with certified rigorous precision bounds.


## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Requirements
Sagee 9.0 is recommended. Furthermore, this project relies on the following packages:

- [Ore Algebra](https://github.com/mkauers/ore_algebra). For better performance, I recommend using the [fastalgexp branch in Marc Mezzarobba's fork](https://github.com/mezzarobba/ore_algebra/tree/fastalgexp).
- [numperiods](https://gitlab.inria.fr/lairez/numperiods) TODO: push coordinates function to a branch of numperiods -- need permissions.
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

## Contact
For any questions, bugs, remarks, please contact [eric.pichon@polytechnique.edu](mailto:eric.pichon@polytechnique.edu)

## Roadmap
Future goals include:
- [ ] Tackling higher dimensional varieties (most notably cubics in P^5).
- [ ] Computing periods of singular varieties.
- [ ] Computing periods of complete intersections.
- [ ] Computing periods of weighed projective hypersurfaces, notably double covers of P^2 ramified along a cubic.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
