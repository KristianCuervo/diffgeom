import sympy as sp
from .manifold import RiemannianManifold
from.curve import Curve

class Field:
    def __init__(self, vector_field: sp.Matrix, space: RiemannianManifold):
        """
        vector_field: a matrix of vectors
        coord_vars: a matrix of coordinates
        """
        self.expr = vector_field
        self.space = space

        self.gradient = self._compute_gradient()

    def _compute_gradient(self):
        """
        Compute the gradient of the vector field.
        """
        grad = []
        for i in range(self.space.dim):
            grad_i = self.expr.diff(self.space.vars[i])
            grad.append(grad_i)
        return sp.Matrix(grad)
    

    def is_killing(self, result=False):
        V = self.expr
        G = self.space.metric
        X = self.space.vars
        K = []
        for i in range(self.space.dim):
            for j in range(self.space.dim):
                K_i = 0
                for k in range(self.space.dim):
                    K_i += (  V[k]*G[i,j].diff(X[k]) 
                            + G[j,k]*V[k].diff(X[i])
                            + G[i,k]*V[k].diff(X[j]))
                    K.append(K_i)
        
        if result == True:
            return K
        elif sum(K) == 0:
            print("Vector is a killing field in G")
            return True
        else:
            print("Vector is not a killing field in G")
            return False
    
    def lie_derivative(self):
            V = self.expr
            G = self.space.metric
            X = self.space.vars
            K = sp.zeros(2)
            for i in range(self.space.dim):
                for j in range(self.space.dim):
                    K_ij = 0
                    for k in range(self.space.dim):
                        K_ij += (  V[k]*G[i,j].diff(X[k]) 
                                + G[j,k]*V[k].diff(X[i])
                                + G[i,k]*V[k].diff(X[j])) 
                    K[i,j] = K_ij
            
            return K

    def is_parallel_along_curve(self, curve: Curve):
        return (sum(self.space.covariant_derivative(vector=self.expr,
                                               curve=curve)) == 0)

class JacobiField(Field):
    def __init__(self, f_s:sp.Matrix, parameter:sp.symbols):
        self.expr = f_s
        self.parameter = parameter
        self.prime = self.expr.diff(parameter)
        self.second_prime = self.prime.diff(parameter)

    def find_sectional_curvature(self):
        if self.expr.subs(self.parameter, 0) == 0:
            print("f(0) = 0 holds")
            condition_1 = True
        if self.prime.subs(self.parameter, 0) == 1:
            print("f'(0) = 1 holds")
            condition_2 = True
        
        if condition_1 and condition_2:
            s = sp.symbols('s')
            sectional_curvature = sp.Function('k')(s)
            eq = sp.Eq(self.second_prime + sectional_curvature*self.expr, 0)
            sol = sp.dsolve(eq, sectional_curvature)
            print(sol)
            return sol
        else:
            if not condition_1:
                print("Condition 1 does not hold:")
            if not condition_2:
                print("Condition 2 does not hold")
            return False
            

        

