import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2025 Problem 3

x1, x2, x3 = sp.symbols('x1 x2 x3')
euclid3d = dg.Space(metric=sp.eye(3),
                    coord_vars=sp.Matrix([x1, x2, x3]))

parameterisation = sp.Matrix([x1*sp.cos(x2),
                              x1*sp.sin(x2),
                              x1])
cone = dg.Surface(embedded_expr=parameterisation,
                  surface_vars=sp.Matrix([x1, x2]),
                  ambient_space=euclid3d)

# 3.1 Show that in x1, x2, the metric is same as for problem 2
print(cone.metric)

# 3.2 Show non-zero christoffel symbols
print(cone.christoffels)
chris122 = cone.christoffels[1,1,0]
chris212 = cone.christoffels[0, 1, 1]
chris221 = cone.christoffels[1, 0, 1]
print(chris122, chris212, chris221)

# 3.3 Constant sectional curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
K = cone.sectional_curvature(X=e1, Y=e2)
print(K) # 0


# 3.4, 3.5 and 3.6 are solved in the pdf for clarity of hand-calculations.
