import sympy as sp
from .space import Space

class Map:
    """
    A map φ from one manifold (domain_space) to another (codomain_space).

    Example:
      domain_space.dim = n
      codomain_space.dim = m

      expr is an m-dimensional sympy Matrix function of n variables,
      i.e. expr : R^n -> R^m

    If n == m and det(Jacobian) != 0, then it can be a diffeomorphism.
    """
    def __init__(self, 
                 expr:sp.Matrix, 
                 domain_space:Space,
                 codomain_space:Space):
        
        # (codomain_vars) = expr(domain_vars)
        self.expr = expr 
        self.domain_space = domain_space
        self.codomain_space = codomain_space

        self.jacobian_matrix = expr.jacobian(domain_space.vars) 

        self.is_diffeomorphism = False
        if self.domain_space.dim == self.codomain_space.dim:
            self.is_diffeomorphism = self._check_diffeomorphism()
    
    def _check_diffeomorphism(self):
        try:
            det_j = self.jacobian_matrix.det()
            return (det_j != 0)
        except (sp.SympifyError, sp.matrices.exceptions.NonSquareMatrixError, TypeError, ValueError):
            return False
    
    def get_inverse_map(self):
        """
        If this is a diffeomorphism (same dim & invertible Jacobian),
        solve codomain_vars = φ(domain_vars) for domain_vars in terms of codomain_vars.
        Return a new Map from codomain_space -> domain_space.
        """
        if not self.is_diffeomorphism:
            #raise ValueError("Map is not a diffeomorphism; cannot compute inverse.")
            pass
        
        # y_i = expr_i(x). Solve for x in terms of y
        eqs = [
            sp.Eq(self.codomain_space.vars[i], self.expr[i])
            for i in range(self.codomain_space.dim)
        ]

        sol = sp.solve(eqs, list(self.domain_space.vars), dict=True)
        if not sol:
            raise ValueError("Could not symbolically solve for inverse.")

        # Take the first solution if multiple
        inv_expr_list = [sol[0][v] for v in self.domain_space.vars]
        inv_expr = sp.Matrix(inv_expr_list)

        # Inverse map goes from codomain_space -> domain_space
        return Map(
            expr=inv_expr,
            domain_space=self.codomain_space,
            codomain_space=self.domain_space
        )

    def pushforward_metric(self):
        """
        Define a new metric on the codomain by pushing the domain
        metric forward under φ.
        1) Check we have a domain metric
        2) Invert the map φ (get φ⁻¹)
        3) Sub x=φ⁻¹(y) into domain_space.metric
        4) Multiply by (Dφ⁻¹)^T (...) (Dφ⁻¹)
        5) Create a new Space object with that metric
        """

        if self.domain_space.metric is None:
            raise ValueError("Domain has no metric to push forward.")

        if not self.is_diffeomorphism:
            raise ValueError("Need diffeo to define a pushforward metric.")

        # 1) Invert φ to get x = φ⁻¹(y)
        inverse_map = self.get_inverse_map()
        J_inv = inverse_map.jacobian_matrix  # D(φ⁻¹)(y

        # 2) Substitutions: x -> φ⁻¹(y)
        x_exprs = inverse_map.expr  # A matrix [x1(y), x2(y), ...]
        sub_dict = {}
        for i in range(self.domain_space.dim):
            sub_dict[self.domain_space.vars[i]] = x_exprs[i]

        domain_metric_in_y = self.domain_space.metric.subs(sub_dict)

        # 3) G_V(y) = J_inv^T * domain_metric_in_y * J_inv
        new_metric = J_inv.T * domain_metric_in_y * J_inv

        # 4) Construct a new Space in codomain coords
        new_space = Space(
            coord_vars= self.codomain_space.vars,
            metric = new_metric
        )

        return new_space