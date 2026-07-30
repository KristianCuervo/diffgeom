# diffgeom

Symbolic Riemannian geometry in Python, built on [SymPy](https://www.sympy.org).

Define a metric — or embed a surface and let the metric be induced — and the
library computes the objects that follow from it exactly, in closed form:
Christoffel symbols, the Riemann curvature tensor, sectional, Ricci and scalar
curvature, covariant derivatives along curves, geodesics, and Killing fields.
Nothing is evaluated numerically; every result comes back as a SymPy expression
you can simplify, solve or differentiate further.

The API deliberately mirrors the notation of the theory, so a written exercise
translates into code more or less line by line.

## Install

```bash
git clone https://github.com/KristianCuervo/differentialGeometry.git
cd differentialGeometry
pip install -e .
```

Requires Python 3.10+, SymPy and NumPy.

## Example

The unit sphere as a surface embedded in Euclidean R³ — its metric is inherited
from the ambient space rather than written down by hand:

```python
import sympy as sp
import diffgeom as dg

x1, x2, x3 = sp.symbols('x1 x2 x3')

euclid3 = dg.Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))

sphere = dg.Surface(
    embedded_expr=sp.Matrix([sp.cos(x2) * sp.cos(x1),
                             sp.cos(x2) * sp.sin(x1),
                             sp.sin(x2)]),
    surface_vars=sp.Matrix([x1, x2]),
    ambient_space=euclid3,
)

e1, e2 = sp.Matrix([1, 0]), sp.Matrix([0, 1])

print(sphere.metric)                        # Matrix([[cos(x2)**2, 0], [0, 1]])
print(sphere.sectional_curvature(e1, e2))   # 1
print(sphere.scalar_curvature())            # 2
```

## What it computes

| Concept | Where |
| --- | --- |
| Metric tensor g(V, W) | `RiemannianManifold.metric_tensor(V, W)` |
| Indicatrix at a point | `RiemannianManifold.find_indicatrix(point)` |
| Christoffel symbols | `RiemannianManifold.christoffels[i, j, m]` |
| Riemann curvature tensor | `RiemannianManifold.riemanns[i, j, k, m]` |
| Curvature operator R(X, Y)Z | `RiemannianManifold.curvature_operator(X, Y, Z)` |
| Curvature tensor R(X, Y, Z, U) | `RiemannianManifold.curvature_tensor(X, Y, Z, U)` |
| Sectional curvature K(X, Y) | `RiemannianManifold.sectional_curvature(X, Y)` |
| Ricci tensor and curvature | `RiemannianManifold.ricci_tensor()`, `.ricci_curvature(X, Y)` |
| Scalar curvature | `RiemannianManifold.scalar_curvature()` |
| Induced metric on a surface | `Surface(...).metric` |
| g-speed and g-length of a curve | `Curve.speed`, `Curve.length(interval, manifold)` |
| Acceleration along a curve | `Curve.acceleration_vector()` |
| Geodesic test | `Curve.is_geodesic()` |
| Gradient of a field | `Field.gradient` |
| Covariant derivative along a curve | `Field.covariant_derivative(curve)` |
| Parallel transport test | `Field.is_parallel_along_curve(curve)` |
| Killing field test | `Field.is_killing()` |
| Lie derivative of the metric | `Field.lie_derivative()` |
| Jacobi fields | `JacobiField.find_sectional_curvature()` |
| Jacobian and diffeomorphism test | `Map.jacobian_matrix`, `Map.is_diffeomorphism` |
| Inverse map and pushforward metric | `Map.get_inverse_map()`, `Map.pushforward_metric()` |

### Classes

- **`RiemannianManifold(vars, metric)`** — the base class. Christoffel symbols
  and the Riemann tensor are computed once on construction and cached.
- **`Space(metric, coord_vars)`** — a manifold given directly by its metric.
- **`EuclideanSpace(n)`** — flat Rⁿ with the identity metric.
- **`Surface(embedded_expr, surface_vars, ambient_space)`** — a surface embedded
  in an ambient manifold; its metric is induced from the embedding.
- **`Curve(curve_expr, parameter, manifold)`** — a parameterised curve.
- **`Field(vector_field, space)`** and **`JacobiField(f_s, parameter)`** — vector
  fields on a manifold.
- **`Map(expr, domain_space, codomain_space)`** — a map between manifolds.

## Layout

```
src/diffgeom/     the library
examples/
  textbook/       worked chapter exercises
  exams/          past exam papers, by year
notebooks/        walkthroughs
tests/            runs every example end to end
```

Every script under `examples/` is a self-contained solution to a real problem,
with the reasoning written out in comments alongside the computation. They run
as-is once the package is installed:

```bash
python examples/exams/2024/p1.py
```

## Tests

```bash
pip install -e ".[dev]"
pytest
```

The suite executes every example and asserts it completes without error, and
checks the curvature routines against metrics with known closed-form answers
(the round sphere, the hyperbolic plane, flat space).

## Origin

Written alongside DTU course 01238 (Differential Geometry). The examples are the
problem sets and past exam papers worked through during the course.

## License

MIT — see [LICENSE](LICENSE).
