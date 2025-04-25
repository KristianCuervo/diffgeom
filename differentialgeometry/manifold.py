import sympy as sp
import numpy as np
from .curve import Curve


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
            return (V.T * self.metric * W)[0, 0]
        else: 
            return (V.T * metric * W)[0, 0]
    
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
                        term += 0.5 * (
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
        return R # Fixes sign convention
    
    def covariant_derivative(self, vector:sp.Matrix, curve: Curve):
        """
        Compute the covariant derivative of a restricted vector
        field V(t) along a curve gamma(t) in the manifold.
        """
        covar = []
        for k in range(self.dim):
            term = vector[k].diff(curve.parameter)
            for i in range(self.dim):
                for j in range(self.dim):
                    term += (vector[j] * 
                             curve.derivative[i] * 
                             self.christoffels[i, j, k].subs(
                                 {self.vars[i]: curve.expr[i] for i in range(self.dim)})) 
            covar.append(term)
        return sp.Matrix(covar)

    def acceleration_vector(self, curve: Curve):
        """
        Compute the acceleration vector of a curve in the manifold.
        """
        return self.covariant_derivative(curve.derivative, curve)
    
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
        K = (R / area)
        
        return K
    
    def scalar_curvature(self):
        """
        Compute the scalar curvature of the manifold.
        """
        S = 0
        for j in range(self.dim):
            for k in range(self.dim):
                for i in range(self.dim):
                    for m in range(self.dim):
                        S += self.riemanns[i, j, k, m] * self.metric_inv[i, m] * self.metric_inv[j, k]
        return S  
    
    def ricci_curvature(self, X:sp.Matrix, Y:sp.Matrix):
        """
        Compute the Ricci curvature of the manifold.
        """
        ric = 0
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for m in range(self.dim):
                        ric += X[i] * Y[j] * self.riemanns[i, j, k, m] * self.metric_inv[k, m]
        return ric
    
    def ricci_tensor(self):
        """
        Compute the Ricci tensor of the manifold.
        """
        return None
    
