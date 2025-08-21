import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2023 Problem 3

x1, x2 = sp.symbols('x1 x2')

space = dg.Space(metric=sp.Matrix([[2*x2**2, 0],
                                   [0, 2+sp.cos(x2)]]),
                coord_vars=[x1, x2])

t = sp.symbols('t')
i_interval = (0, 1)

curve = dg.Curve(curve_expr=sp.Matrix([t, t]),
                 parameter=t,
                 manifold=space)

curve_speed = curve.speed
print(curve_speed)

length = curve.length(interval=i_interval, manifold=space)
print(length) # Answered as an integral