import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Problem 7.45

x1, x2 = sp.symbols('x1 x2', real=True)
x3 = sp.symbols('x3', real=True, positive=True)
f = (1/(x3))

G = f*f*sp.eye(3)

space = dg.Space(metric=G, coord_vars=sp.Matrix([x1, x2, x3]))

K = space.sectional_curvature(X=sp.Matrix([1,0,0]), Y=sp.Matrix([0,1,0]))
print(K)

