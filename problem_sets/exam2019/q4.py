import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2019 Problem 4
x1, x2 = sp.symbols('x1 x2')
space = dg.Space(metric=(x1**2)*sp.eye(2),
                 coord_vars=sp.Matrix([x1, x2]))

e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
R1221 = space.riemanns[0, 1, 1, 0]
print(R1221)

# find sectioanl curvatureK
K = space.sectional_curvature(X=e1, Y=e2)
print(K)
