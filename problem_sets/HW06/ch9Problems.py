import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# 9.19

# AMBIENT SPACE SETUP
y1, y2 = sp.symbols('y1, y2')
y3 = sp.symbols('y3', positive=True)

ambient = dg.Space(metric=((1/y3)**2)*sp.eye(3),
                         coord_vars=[y1, y2, y3])

# SURFACE INSIDE OF AMBIENT SPACE
x1, x2 = sp.symbols('x1, x2')
theta = sp.symbols('theta')
phi = sp.Matrix([x1, x2*sp.cos(theta), 1+x2*sp.sin(theta)])
surface_phi = dg.Surface(embedded_expr=phi, surface_vars=[x1, x2],
                         ambient_space=ambient)

# Show induced metric g is of certain form
print(surface_phi.metric) # True

# 9.20 Show that (U, g, levi) has constant intrinsic sectional curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])

intrinsic_curvature = surface_phi.sectional_curvature(X=e1, Y=e2)
print(intrinsic_curvature) # True

