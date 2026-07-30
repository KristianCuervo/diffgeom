"""
Curvature checks against metrics whose answers are known in closed form.

These guard the index conventions in ricci_tensor and scalar_curvature.
Both were previously contracted incorrectly: the Ricci contraction vanished
identically for every metric, and the scalar curvature raised an index that
was already raised, giving 1 + 1/sin^2(x) for the unit 2-sphere instead of 2.
"""
import sympy as sp
import pytest

import diffgeom as dg

x1, x2, x3 = sp.symbols('x1 x2 x3')
e1_2d, e2_2d = sp.Matrix([1, 0]), sp.Matrix([0, 1])


def unit_sphere():
    """The round 2-sphere: K = 1, Ric = g, S = 2."""
    return dg.Space(metric=sp.Matrix([[1, 0], [0, sp.sin(x1) ** 2]]),
                    coord_vars=sp.Matrix([x1, x2]))


def hyperbolic_plane():
    """The upper half-plane model: K = -1, Ric = -g, S = -2."""
    return dg.Space(metric=sp.Matrix([[1 / x2 ** 2, 0], [0, 1 / x2 ** 2]]),
                    coord_vars=sp.Matrix([x1, x2]))


def flat_plane():
    return dg.Space(metric=sp.eye(2), coord_vars=sp.Matrix([x1, x2]))


def flat_space():
    return dg.Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))


@pytest.mark.parametrize("build, expected", [
    (unit_sphere, 2),
    (hyperbolic_plane, -2),
    (flat_plane, 0),
    (flat_space, 0),
])
def test_scalar_curvature(build, expected):
    assert sp.simplify(build().scalar_curvature() - expected) == 0


@pytest.mark.parametrize("build, factor", [
    (unit_sphere, 1),        # Ric = +g
    (hyperbolic_plane, -1),  # Ric = -g
    (flat_plane, 0),         # Ric = 0
])
def test_ricci_tensor_is_proportional_to_metric(build, factor):
    space = build()
    assert sp.simplify(space.ricci_tensor() - factor * space.metric) == sp.zeros(2, 2)


@pytest.mark.parametrize("build, expected", [
    (unit_sphere, 1),
    (hyperbolic_plane, -1),
    (flat_plane, 0),
])
def test_sectional_curvature(build, expected):
    K = build().sectional_curvature(e1_2d, e2_2d)
    assert sp.simplify(K - expected) == 0


def test_ricci_curvature_matches_tensor():
    space = unit_sphere()
    ric = space.ricci_tensor()
    for X in (e1_2d, e2_2d):
        for Y in (e1_2d, e2_2d):
            assert sp.simplify(space.ricci_curvature(X, Y) - (X.T * ric * Y)[0]) == 0


def test_christoffel_symbols_stay_exact():
    """Rational(1, 2) rather than 0.5 keeps symbolic output free of floats."""
    space = unit_sphere()
    for i in range(space.dim):
        for j in range(space.dim):
            for m in range(space.dim):
                assert not space.christoffels[i, j, m].atoms(sp.Float)


def test_induced_metric_on_embedded_sphere():
    euclid3 = dg.Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))
    sphere = dg.Surface(
        embedded_expr=sp.Matrix([sp.cos(x2) * sp.cos(x1),
                                 sp.cos(x2) * sp.sin(x1),
                                 sp.sin(x2)]),
        surface_vars=sp.Matrix([x1, x2]),
        ambient_space=euclid3,
    )
    expected = sp.Matrix([[sp.cos(x2) ** 2, 0], [0, 1]])
    assert sp.simplify(sphere.metric - expected) == sp.zeros(2, 2)
    assert sp.simplify(sphere.sectional_curvature(e1_2d, e2_2d) - 1) == 0


def test_meridian_is_a_geodesic():
    # x1 is the polar angle, so the curve (t, 0) is a meridian: a great
    # circle, and therefore a geodesic.
    t = sp.symbols('t')
    meridian = dg.Curve(curve_expr=sp.Matrix([t, 0]),
                        parameter=t,
                        manifold=unit_sphere())
    assert meridian.is_geodesic()


def test_euclidean_space_is_flat():
    space = dg.EuclideanSpace(3)
    assert space.scalar_curvature() == 0
    assert space.ricci_tensor() == sp.zeros(3, 3)
