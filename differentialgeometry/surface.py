import sympy as sp
from .manifold import RiemannianManifold

class Surface(RiemannianManifold):
    def __init__(self, embedded_expr: sp.Matrix, surface_vars=sp.Matrix):
        """
        embedded_expr: r(u,v) in R^3
        surface_vars: (u, v)
        """
        Jr = embedded_expr.jacobian(surface_vars)
        induced_metric = Jr.T * Jr

        super().__init__(vars=surface_vars, metric=induced_metric)
        self.expr = embedded_expr