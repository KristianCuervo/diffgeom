import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2019 Problem 1

x1, x2 = sp.symbols('x1, x2')
G = sp.Matrix([[1 + x1**2, 1],
               [1, 2 + x2**2]])
space = dg.Space(metric=G, coord_vars=sp.Matrix([x1, x2]))

p = (0, 0)
t = sp.symbols('t')
E = sp.Matrix([sp.cos(t), 2*sp.sin(t)])

# 1.1 Show that E(t) is NOT the indicatrix
""""
For I(p), we require that E * Gu * E.T = 1 at p
"""
indicatrix_test = space.metric_tensor(V=E, W=E).subs({x1: p[0],
                                                     x2: p[1]})
print(sp.trigsimp(indicatrix_test))
"""
g(E, E) = 2*sin(2*t) - 7*cos(2*t)/2 + 9/2
Which is not equal to one for all t on [-pi, pi]
"""

# 1.2 Find indicatrix
i, it, t = space.find_indicatrix(point=p)
print((sp.trigsimp(it)))

# 1.3 Find g-length of the curve eta
eta = dg.Curve(curve_expr=sp.Matrix([t, 1]),
               parameter=t,
               manifold=space)
interval = (-7, 7)
g_length = eta.length(interval=interval, manifold=space)
print("G-length is equal to:", g_length)

# 1.4 Show that eta is not a geodesic in M2
print(eta.is_geodesic(space=space))

# 1.5 Show that eta is not the g-shortest curve between its endpoints


