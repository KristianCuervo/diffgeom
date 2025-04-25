import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp


# Problem 7.47
x1, x2 = sp.symbols('x1 x2', real=True)
x3 = sp.symbols('x3', real=True, positive=True)
f = (1/(x3))

G = sp.eye(3)
G[2,2] = f*f

space = dg.Space(metric=G, coord_vars=sp.Matrix([x1, x2, x3]))

# compute scalar curvature
S = space.scalar_curvature()
print(S)

# Ricci curvature in e1, e2, e3

e1 = sp.Matrix([1, 0, 0])
e2 = sp.Matrix([0, 1, 0])
e3 = sp.Matrix([0, 0, 1])

re1 = space.ricci_curvature(e1, e1)
print(re1)

re2 = space.ricci_curvature(e2, e2)
print(re2)

re3 = space.ricci_curvature(e3, e3)
print(re3)

# Three spaces

ke1e2 = space.sectional_curvature(e1, e2)
print(ke1e2)
ke1e3 = space.sectional_curvature(e1, e3)
print(ke1e3)
ke2e3 = space.sectional_curvature(e2, e3)
print(ke2e3)


