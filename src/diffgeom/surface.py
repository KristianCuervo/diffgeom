import sympy as sp
from .manifold import RiemannianManifold

class Surface(RiemannianManifold):
    def __init__(self, embedded_expr: sp.Matrix, surface_vars=sp.Matrix, ambient_space=RiemannianManifold):
        """
        embedded_expr: r(x1, x2,...,xn) in R^N
        surface_vars: e.g x1, x2,...,xn

        """
        self.expr = embedded_expr
        Jr = self.expr.jacobian(surface_vars)
        
        G_sub = ambient_space.metric.subs({
            ambient_space.vars[i]: self.expr[i]
            for i in range(ambient_space.dim)
        })

        induced_metric = sp.simplify(ambient_space.metric_tensor(V=Jr, W=Jr, metric=G_sub))

        super().__init__(metric=induced_metric, vars=surface_vars)