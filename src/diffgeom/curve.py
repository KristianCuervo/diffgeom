import sympy as sp
from .space import Space

class Curve:
    def __init__(self, curve_expr: sp.Matrix, parameter: sp.Symbol, manifold: Space):
        """
        curve_expr: gamma(t), e.g. a vector in R^n
        parameter: the symbol t
        """
        self.expr = curve_expr
        self.parameter = parameter
        
        self.derivative = self.expr.diff(self.parameter)

        self.manifold = manifold
        self.speed = self._compute_prime(manifold)

    def _compute_prime(self, manifold: Space):
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

    def length(self, interval: tuple, manifold: Space):
        """
        Compute ∫ sqrt( (dγ/dt)^T * G(γ(t)) * (dγ/dt) ) dt 
        from t0 to t1, if manifold has a metric.
        """
        t0, t1 = interval
        return sp.integrate(self.speed, (self.parameter, t0, t1))
    
    
    def is_geodesic(self) -> bool:
        """
        A curve is a geodesic when its acceleration vector vanishes.
        """
        return all(sp.simplify(a) == 0 for a in self.acceleration_vector())

    def acceleration_vector(self):
        """
        Compute the covariant derivative of a restricted vector
        field V(t) along a curve gamma(t) in the manifold.
        """
        V = self.derivative
        gammaprime = self.derivative
        man = self.manifold
        along_curve = {man.vars[n]: self.expr[n] for n in range(man.dim)}
        covar = []
        for k in range(man.dim):
            covar_i = 0
            covar_i += V.diff(self.parameter)[k]
            for i in range(man.dim):
                for j in range(man.dim):
                    covar_i += (V[i] *
                             gammaprime[j] *
                             man.christoffels[i, j, k].subs(along_curve))
            covar.append(sp.trigsimp(covar_i))
        return sp.Matrix(covar)
