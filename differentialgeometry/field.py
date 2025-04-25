import sympy as sp
from .space import Space
from.curve import Curve

class Field:
    def __init__(self, vector_field: sp.Matrix, space: Space):
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
    
