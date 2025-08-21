import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2022 Problem 1

x1, x2 = sp.symbols('x1 x2')
G = sp.Matrix([[2 + x1**2, 0],
               [0 , 1]])

space = dg.Space(metric=G,
                 coord_vars=sp.Matrix([x1, x2]))

# 1.1 
"""
The condition must be met that 
I(p)=1= V * Gv * V.T 
"""
v1 = sp.Matrix([1, 0])
v2 = sp.Matrix([0, 1])
v3 = sp.Matrix([1/sp.sqrt(2), 0])
v4 = sp.Matrix([0, 1/sp.sqrt(2)])
v_list = [v1, v2, v3, v4]
for i, v_i in enumerate(v_list):
    print(f"I_{i+1} = ", space.metric_tensor(V=v_i, W=v_i,
                                             metric=space.metric.subs({x1: 0,
                                                                       x2: 0})))
# This is only equal to 1 for all x for v2
# FALSE: It is actually important to check these at (0, 0)
# Then it is v2, v3 which are equal to 1

# 1.2 Find all non-zero christoffel symbols
chris = space.christoffels
#print(chris)
"""
The only non-zero christoffel symbol as seen from the print is 
chris_1, 1 = 1.0*x1/(x1**2 + 2)
"""
print(chris[0][0][0]) # This is the only non-zero christoffel symbol

# 1.3  Find the arc-length s(t) along alpha(t) 
t = sp.symbols('t')
alpha = dg.Curve(curve_expr=sp.Matrix([t, 1]),
                 parameter=t, 
                 manifold=space)
s_arclength_function = sp.Integral(alpha.speed)
#print(s_arclength_function) # from -1 to t

# 1.4 Determine whether alpha can be reparameterised as a geodesic
covar_alpha = space.covariant_derivative(vector=alpha.derivative,
                                         curve=alpha)
#print(covar_alpha) # This must be equal to zero. But there is a non-zero term here.
"""
A pre-geodesic is found by
covar_alpha = rho(t) * alpha.derivaitve(t)
Thus we try to find rho(t)
"""
rho = covar_alpha[0]/alpha.derivative[0]
"""
rho can be found and thus alpha is pre-geodesic.
Thus alpha can be reparametersed as a geodesic.
"""

# 1,5 Find tangent vector and acceleration of beta
beta = dg.Curve(curve_expr=sp.Matrix([t, t**2]),
                parameter=t,
                manifold=space)

tangent_vector = beta.derivative.subs(t, 0)
print(tangent_vector)

acceleration_beta = space.acceleration_vector(curve=beta).subs(t, 0)
print(acceleration_beta)

# Both alpha and beta are paths in U, which has shorter g-length between end-points
point_p = (-1, 1) # t = -1
point_q = (1, 1) # t = 1

t_interval = (-1, 1)

glength_alpha = alpha.length(interval=t_interval,
                             manifold=space)
glength_beta = beta.length(interval=t_interval,
                           manifold=space)

print(glength_alpha)
print(glength_beta)
if glength_alpha > glength_beta:
    print("beta is shorter")
else:
    print("alpha is shorter")

e1 = dg.Field(vector_field=sp.Matrix([1, 0]),
              space=space)
e2 = dg.Field(vector_field=sp.Matrix([0, 1]),
              space=space)

print(e1.is_killing()) # e1 is not a killing field in G
print(e2.is_killing()) # e2 IS a killing field in G

# 1.8 Explain why 1.7 is expected
"""
The g-metric is only dependent on x1, which means that the length
of e1 is not conserved along a path in the space U.

In contrast, e2 follows in the direction of e2. The space U
does not vary its metric along any straight path in the direction
of e2. Thus the length of the vectors in the vectorfield are conserved
along their path.
"""