import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2022 Problem 4

x1, x2 = sp.symbols('x1, x2')
y1, y2 = sp.symbols('y1, y2')

space_U = dg.Space(metric=sp.eye(2),
                   coord_vars=sp.Matrix([x1, x2]))

space_V = dg.Space(metric=sp.eye(2),
                   coord_vars=sp.Matrix([y1, y2]))

phi = dg.Map(expr=sp.Matrix([x1**3, x2**3]),
             domain_space=space_U,
             codomain_space=space_V)

# 4.1 Show that phi: U->V is a bijection
"""
Phi is a bijection if the inverse is well defined.
Furthermore, the partial derivatives exist and the
jacobian matrix and the jacobian of the inverse map is well-defined.
"""
phi_inv = phi.get_inverse_map()
print(phi_inv.expr) # well-defined inverse expression

print(phi.jacobian_matrix) # well-defined and exists on x1, x2
#print(phi_inv.jacobian_matrix) # well-defined and exists on y1, y2 where y1 or y2 != 0

# 4.2 Find the jacobian of the inverse map
"""
This is yet again given by:
"""
print(phi_inv.jacobian_matrix)

# 4.3 Explain why U-->V is NOT a diffeomorphism
"""
For phi to be a diffeomorphism, then its jacobian and 
the jacobian of the inverse function must be defined for all 
values in the space.

However, when looking at the jacobian of the inverse map, it is clear
that it is not defined on the axis of the space:
when y1=0 or y2 = 0.
"""

# 4,4 Explain why phi+ is a diffeomorphism between U+ and V+.
"""
U+ and V+ only take into account the first quadrant of the respective spaces,
when x1, x2 > 0 and y1, y2 > 0. For these values, then the argument in 
4.3 does not hold the same way. Now the inverse map is well-defined for
all values on the entire domains and the points y1=0 and y2=0 are not in the spaces
U+ and V+
"""

# 4.5 Get the metric matrix function as a function of y1, y2
space_V_new = phi.fromUtoV()
print(space_V_new.metric)

# 4.6 find the sectional curvature in V
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0 ,1])
print(space_V_new.sectional_curvature(X=e1, Y=e2)) # 0

# 4.7 find g-length of path 
t = sp.symbols('t')
interval = (1, 8)
gamma = dg.Curve(curve_expr=sp.Matrix([t, 1]),
                 parameter=t,
                 manifold=space_V_new)

glength_gamma = gamma.length(interval=interval,
                             manifold=space_V_new)
print(glength_gamma) # 1



