from __future__ import annotations
import sympy as sp

class Curve:
    def __init__(self, curve_expr: sp.Matrix, parameter: sp.Symbol):
        """
        curve_expr: gamma(t), e.g. a vector in R^n
        parameter: the symbol t
        """
        self.expr = curve_expr
        self.parameter = parameter

        self.derivative = self._compute_prime()

    def _compute_prime(self):
        """
        Compute the derivative of the curve with respect to the parameter.
        """
        return self.expr.diff(self.parameter)

    def length(self, interval: tuple, manifold: "RiemannianManifold"):
        """
        Compute ∫ sqrt( (dγ/dt)^T * G(γ(t)) * (dγ/dt) ) dt 
        from t0 to t1, if manifold has a metric.
        """
        t0, t1 = interval
        dgamma_dt = self.curve_expr.diff(self.parameter)
        
        subs_dict = {
            manifold.vars[i]: self.expr[i]
            for i in range(manifold.dim)
        }
        G_sub = manifold.metric.subs(subs_dict)
        
        integrand = sp.sqrt((dgamma_dt.T * G_sub * dgamma_dt)[0, 0])
        return sp.integrate(integrand, (self.parameter, t0, t1))
    
    
    def is_geodesic(self, space: "RiemannianManifold") -> bool:
        acceleration = space.acceleration_vector(curve=self)
        return all(a == 0 for a in acceleration)
    