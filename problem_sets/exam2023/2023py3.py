import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# 2022 Problem 3

x1, x2 = sp.symbols('x1, x2')
G = sp.Matrix([[2 + x2**2, 0],
               [0, 2 + sp.cos(x2)]])

space = dg.Space(metric=G,
                 coord_vars=sp.Matrix([x1, x2]))

t = sp.symbols('t')
alpha = dg.Curve(curve_expr=sp.Matrix([t, t]),
                 parameter=t,
                 manifold=space)

# 3.1 Find the g-speed of the curve as a function of t

print(alpha.speed)

# 3.2 g-length from t=0 to t=1
t_interval = (0, 1)
g_length = alpha.length(interval=t_interval,
                        manifold=space)
print(g_length)

# 3.3 Find sectional curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
k_sectional = sp.simplify(sp.trigsimp(space.sectional_curvature(X=e1, Y=e2)))
print(k_sectional)

# 3.4 geodesic cirlce of radiu e. 
"""
We take a look at the behaviour of the metric around the origin, and investigate
whether the metric will stretch or condense the space.
This stretching or condensing is relative to the value of the eigenvalues, 
which represent the 
"""
G_origin = G.subs({x1:0,
                   x2:0})
print(G_origin)
# G = [[2, 0], [0, 3]]
"""
We test this using eucli2d as a space and also
an origin space which is the metric evaluated at the origin.
"""
origin_space = dg.Space(metric=sp.Matrix([[2, 0],
                                         [0, 3]]),
                        coord_vars=sp.Matrix([x1, x2]))

euclid2d = dg.Space(metric=sp.eye(2),
                    coord_vars=sp.Matrix([x1, x2]))
circle = dg.Curve(curve_expr=sp.Matrix([sp.cos(t), sp.sin(t)]),
                  parameter=t,
                  manifold=euclid2d)
t_range = (-sp.pi, sp.pi)
euclidean_perimeter = circle.length(interval=t_range,
                                    manifold=euclid2d)

circle_stretched =  dg.Curve(curve_expr=sp.Matrix([sp.cos(t), sp.sin(t)]),
                  parameter=t,
                  manifold=origin_space)
stretched_perimeter = circle_stretched.length(interval=t_range,
                                    manifold=origin_space)
print(euclidean_perimeter)
print(stretched_perimeter.evalf(5)) # MORE THAN 2PI???

# Alternative, we can look at the curvature. 
k_origin = k_sectional.subs({x1:0,
                             x2:0})
print(sp.nsimplify(k_origin)) # -1/6 
"""
The curvature is negative and thus the length of the 
small sectional geodesic curvatures will be larger than 2pi*epsilon.
This is seen as a result from theorem 6.29.
"""
