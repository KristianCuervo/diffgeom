# diffgeom

Symbolic Riemannian geometry in Python, built on [SymPy](https://www.sympy.org).

You give it a metric, or you embed a surface and let it work out the induced
metric. From there it computes Christoffel symbols, the Riemann curvature
tensor, sectional, Ricci and scalar curvature, covariant derivatives along
curves, geodesics and Killing fields. Everything stays symbolic. Results come
back as SymPy expressions, so you can keep simplifying, solving or
differentiating them.

The class and method names follow the notation used in the theory, which means
a written exercise usually translates into code close to line by line.

## Install

```bash
git clone https://github.com/KristianCuervo/diffgeom.git
cd diffgeom
pip install -e .
```

Requires Python 3.10+, SymPy and NumPy.

## Example

The unit sphere, treated as a surface embedded in Euclidean R³. Its metric is
inherited from the ambient space rather than written down by hand:

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
| Metric tensor g(V, W) | `Space.metric_tensor(V, W)` |
| Indicatrix at a point | `Space.find_indicatrix(point)` |
| Christoffel symbols | `Space.christoffels[i, j, m]` |
| Riemann curvature tensor | `Space.riemanns[i, j, k, m]` |
| Curvature operator R(X, Y)Z | `Space.curvature_operator(X, Y, Z)` |
| Curvature tensor R(X, Y, Z, U) | `Space.curvature_tensor(X, Y, Z, U)` |
| Sectional curvature K(X, Y) | `Space.sectional_curvature(X, Y)` |
| Ricci tensor and curvature | `Space.ricci_tensor()`, `.ricci_curvature(X, Y)` |
| Scalar curvature | `Space.scalar_curvature()` |
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

**`Space(metric, coord_vars)`** is a Riemannian manifold given directly by its
metric, and the base class for the rest. It computes the Christoffel symbols and
the Riemann tensor once on construction and caches them.
**`EuclideanSpace(n)`** is flat Rⁿ with the identity metric.

**`Surface(embedded_expr, surface_vars, ambient_space)`** is a surface embedded
in an ambient manifold. Its metric is induced from the embedding.

**`Curve(curve_expr, parameter, manifold)`** is a parameterised curve.

**`Field(vector_field, space)`** and **`JacobiField(f_s, parameter)`** are vector
fields on a manifold.

**`Map(expr, domain_space, codomain_space)`** is a map between manifolds.

## Layout

```
src/diffgeom/     the library
examples/
  textbook/       worked chapter exercises
  exams/          past exam papers, by year
```

Each script under `examples/` is a self-contained solution to a real problem,
with the reasoning written out in comments next to the computation. They run
as-is once the package is installed:

```bash
python examples/exams/2024/p1.py
```

## Origin

I wrote this while taking DTU course 01238 (Differential Geometry), taught by
[Steen Markvorsen](http://www2.mat.dtu.dk/people/S.Markvorsen/). The notation and
the structure of the classes follow the course material. The examples are the
problem sets and past exam papers I worked through during the course.

## License

MIT, see [LICENSE](LICENSE).
