from __future__ import annotations
import sympy as sp

class Curve:
    def __init__(self, curve_expr: sp.Matrix, parameter: sp.Symbol, manifold: "RiemannianManifold"):
        """
        curve_expr: gamma(t), e.g. a vector in R^n
        parameter: the symbol t
        """
        self.expr = curve_expr
        self.parameter = parameter

        self.derivative = self.expr.diff(self.parameter)
        self.speed = self._compute_prime(manifold)

    def _compute_prime(self, manifold: "RiemannianManifold"):
        """
        Compute the derivative of the curve with respect to the parameter and g-metric.
        """
        dgamma_dt = self.expr.diff(self.parameter)
        
        subs_dict = {
            manifold.vars[i]: self.expr[i]
            for i in range(manifold.dim)
        }
        G_sub = manifold.metric.subs(subs_dict)
        
        g_speed = sp.sqrt(manifold.metric_tensor(V=dgamma_dt, W=dgamma_dt, metric=G_sub)[0])

        return g_speed

    def length(self, interval: tuple, manifold: "RiemannianManifold"):
        """
        Compute ∫ sqrt( (dγ/dt)^T * G(γ(t)) * (dγ/dt) ) dt 
        from t0 to t1, if manifold has a metric.
        """
        t0, t1 = interval
        return sp.integrate(self.speed, (self.parameter, t0, t1))
    
    
    def is_geodesic(self, space: "RiemannianManifold") -> bool:
        acceleration = space.acceleration_vector(curve=self)
        return all(a == 0 for a in acceleration)
    
    