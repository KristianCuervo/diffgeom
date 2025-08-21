import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2025 Problem 1
y1, y2, y3 = sp.symbols('y1 y2 y3')
H_v = (1/((y1**2 + y2**2 + y3**2)**2))*sp.eye(3)

space_m = dg.Space(metric=H_v,
                   coord_vars=sp.Matrix([y1, y2, y3]))


t = sp.symbols('t', real=True)
gamma = dg.Curve(curve_expr=sp.Matrix([t, 0, 1]),
                 parameter=t,
                 manifold=space_m)

# 1.1 Find the length of the curve segment gamma
t_interval = (0, 1)
length = gamma.length(interval=t_interval,
                      manifold=space_m)
print("length of gamma = ", length) # pi/4

# 1.2 Find the three christoffel symbols along gamma as functions of the parameter t
# chris[k, 11](t, 0, 1): k = 1,2,3
chris111_gamma = space_m.christoffels[0,0,0].subs({space_m.vars[i]: gamma.expr[i] for i in range(space_m.dim)})
chris211_gamma = space_m.christoffels[0,0,1].subs({space_m.vars[i]: gamma.expr[i] for i in range(space_m.dim)})
chris311_gamma = space_m.christoffels[0,0,2].subs({space_m.vars[i]: gamma.expr[i] for i in range(space_m.dim)})

print("chris111 = ",chris111_gamma)
print("chris211 = ",chris211_gamma)
print("chris311 = ",chris311_gamma)

# 1.3 Show that gamma is not a pre-geodesic curve in M
acc_gamma = space_m.acceleration_vector(curve=gamma)
print(acc_gamma)
print(gamma.derivative)
"""
Gamma is not a pre-geodesic because it is not possible to solve
for the definition 3.46 for a geodesic and find a smooth
function rho of t. 

Comparing acc_gamma and gamma.derivative, we see that 
the third entry of gamma'(t) is equal to zero, and thus it can't
be multiplied by a rho(t) to become the third entry of 
the acceleration vector of gamma.
"""

# 1.4 Find Gu
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
""""
The parameterisation N can be seen as a surface which maps
from U=R2 to the ambient space V previously constructed.

"""
surface_n = dg.Surface(embedded_expr=sp.Matrix([x1, x2, 1]),
                       surface_vars=sp.Matrix([x1, x2]),
                       ambient_space=space_m)
Gu = surface_n.metric
print(Gu) 

# 1.5 Find the sectional curvature Kp of (U, gv) in terms of x1 x2
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])

K_p = surface_n.sectional_curvature(X=e1, Y=e2)
print("Sectional curvature = ", K_p) # 4.0

# 1.6 Show that eta is a pre-geodesic curve in (U, gu)
eta = dg.Curve(curve_expr=sp.Matrix([t, 0]),
               parameter=t,
               manifold=surface_n)
acc_eta = surface_n.acceleration_vector(curve=eta)
print(acc_eta)
print(gamma.derivative)
"""
In contrast to gamma, the acceleration vector can actually
be constructed using a smooth function rho(t) and 
its derivative according to definition 3.47
"""
rho = acc_eta[0] / gamma.derivative[0]
print(rho) 
"""
The function is also clearly smooth because as the denominator
is equal to t**2 + 1, then there is no case where it is equal to
zero and there is a discontinuity.
"""

