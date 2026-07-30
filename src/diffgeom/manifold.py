import sympy as sp
import numpy as np


def _simplify(expr):
    """
    Simplify a curvature expression. Contractions routinely produce mixed
    tan/sin/cos forms that sp.simplify alone leaves untouched, so rewrite
    everything in terms of sin first.
    """
    expr = sp.sympify(expr)
    try:
        return sp.simplify(expr.rewrite(sp.sin))
    except (AttributeError, TypeError, ValueError):
        return sp.simplify(expr)


class RiemannianManifold:
    def __init__(self, vars:sp.Matrix, metric:sp.Matrix=None):
        """
        Generic n-dimensional manifold with optional metric.
        Used as a superclass.
        vars = coordinates in manifold
        metric = metric matrix G of manifold (if given).
        """

        self.vars = vars
        self.dim = len(vars)
        self.metric = metric
        self.metric_inv = None
        if metric is None:
            raise ValueError(
                "A metric is required to compute the geometry of the manifold."
            )
        if metric.det() != 0:
            self.metric_inv = self.metric.inv()
        
        # Precompute geometry
        self.christoffels = self._compute_christoffel_symbols()
        self.riemanns = self._compute_riemann_tensor()


    def metric_tensor(self, V:sp.Matrix, W:sp.Matrix, metric=None):
        """
        g(V, W) = V^T * metric * W
        """
        if metric is None:  
            return (V.T * self.metric * W)
        else: 
            return (V.T * metric * W)
        
    def find_indicatrix(self, point):
        G_at_point = self.metric.subs({self.vars[i]: point[i] for i in range(self.dim)})
        V = sp.Matrix([x for x in self.vars])
        Ip = self.metric_tensor(V=V, W=V, metric=G_at_point)
        t = sp.symbols('t')
        polar = sp.Matrix([sp.cos(t), sp.sin(t)])
        Ipt = Ip.subs({self.vars[i]: polar[i] for i in range(self.dim)})
        return Ip, Ipt, t
    
    def _compute_christoffel_symbols(self):
        """
        Returns a sympy 3D array of Γ^m_{ij}.
        Access using christoffel[i,j,m]
        """
        chris = sp.MutableDenseNDimArray.zeros(self.dim, self.dim, self.dim)
        #chris = sp.zeros(self.dim, self.dim, self.dim)
        for i in range(self.dim):
            for j in range(self.dim):
                for m in range(self.dim):
                    term = 0
                    for l in range(self.dim):
                        term += sp.Rational(1, 2) * (
                            self.metric[j, l].diff(self.vars[i]) +
                            self.metric[l, i].diff(self.vars[j]) -
                            self.metric[i, j].diff(self.vars[l])
                        ) * self.metric_inv[m, l]
                    chris[i, j, m] = sp.simplify(term)
        return chris

    def _compute_riemann_tensor(self) -> np.array:
        """
        Returns a 4D numpy array R^m_{ijk}.
        Access using R[i,j,k,m]
        """
        R = np.zeros((self.dim, self.dim, self.dim, self.dim), dtype=object)
        chris = self.christoffels
        X = self.vars
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for m in range(self.dim):
                        term = chris[j, k, m].diff(X[i]) - chris[i, k, m].diff(X[j])
                        for s in range(self.dim):
                            term += chris[j, k, s]*chris[i, s, m] - chris[i, k, s]*chris[j, s, m]
                        R[i, j, k, m] = sp.simplify(term)
        return R
    
    def curvature_operator(self, X:sp.Matrix, Y:sp.Matrix, Z:sp.Matrix):
        """
        Compute the curvature tensor R(X,Y)Z = R^m_{ijk} X^i Y^j Z^k
        """
        R = []
        for m in range(self.dim):
            m_term = 0
            for i in range(self.dim):
                for j in range(self.dim):
                    for k in range(self.dim):
                        m_term += self.riemanns[i, j, k, m] * X[i] * Y[j] * Z[k]
            R.append(m_term)
        return sp.Matrix(R) 

    def curvature_tensor(self, X:sp.Matrix, Y:sp.Matrix, Z:sp.Matrix, U:sp.Matrix):
        return self.metric_tensor(V = self.curvature_operator(X, Y, Z), W = U) 
    
    def sectional_curvature(self, X:sp.Matrix, Y:sp.Matrix, metric:sp.Matrix = None):
        """
        Compute sectional curvature 
        K(X,Y)= R(X,Y,Y,X)/Area2g(X,Y)
        """
        if metric is None or metric == self.metric:
            metric = self.metric
        else:
            metric = metric.subs({self.vars[i]: X[i] for i in range(self.dim)})

        # Compute the area form
        xnorm2 = self.metric_tensor(X, X, metric)
        ynorm2 = self.metric_tensor(Y, Y, metric)
        gXY2 = self.metric_tensor(X, Y, metric)
        area = xnorm2 * ynorm2 - gXY2**2

        R = self.curvature_tensor(X=X, Y=Y, Z=Y, U=X)
        # Compute the sectional curvature
        K = (R[0] / area[0])
        
        return K
    
    def scalar_curvature(self):
        """
        Compute the scalar curvature S = g^{jk} Ric_{kj}.
        """
        ric = self.ricci_tensor()
        S = 0
        for j in range(self.dim):
            for k in range(self.dim):
                S += self.metric_inv[j, k] * ric[k, j]
        return _simplify(S)


    def ricci_curvature(self, X:sp.Matrix, Y:sp.Matrix):
        """
        Compute the Ricci curvature Ric(X, Y) = X^i Y^j Ric_{ij}.
        """
        return (X.T * self.ricci_tensor() * Y)[0]

    def ricci_tensor(self):
        """
        Compute the Ricci tensor Ric_{kj} = R^i_{kij} as a dim x dim
        sympy Matrix, by contracting the upper index of the Riemann
        tensor against its first lower index.

        Note self.riemanns stores R[i, j, k, m] = R^m_{kij}, so the
        contraction is over the first and last stored slots.
        """
        ric = sp.zeros(self.dim, self.dim)
        for k in range(self.dim):
            for j in range(self.dim):
                term = 0
                for i in range(self.dim):
                    term += self.riemanns[i, j, k, i]
                ric[k, j] = _simplify(term)
        return ric


class EuclideanSpace(RiemannianManifold):
    def __init__(self, n:int):
        """
        The flat n-dimensional Euclidean space with coordinates x1..xn
        and the identity metric.
        """
        coord_vars = sp.Matrix(sp.symbols(f'x1:{n + 1}', real=True))
        super().__init__(vars=coord_vars, metric=sp.eye(n))