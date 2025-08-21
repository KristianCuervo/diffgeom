import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from differentialgeometry import Surface, Curve, Space, Map
import sympy as sp

"""
Exercise 1.29
"""
x1, x2 = sp.symbols('x1 x2')
x3 = sp.symbols('x3')
euclid3d = Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))

# surface defined
surfaceR = Surface(embedded_expr=sp.Matrix([x1, x2 , x1**2 + x2**2]), surface_vars=sp.Matrix([x1, x2]))
print(surfaceR.metric) # a) Gu is correct


# Diffeomorphism phi
phiMap = Map(expr=sp.Matrix([sp.sqrt(x1**2 + x2**2), sp.arg(x1 + x2*sp.I)]),
             domain_space=euclid3d, codomain_space=surfaceR)

# Curves
interval = (-sp.pi, +sp.pi)
t = sp.symbols('t')

# Using formula 1.48
gamma = Curve(curve_expr=sp.Matrix([sp.cos(t), sp.sin(t)]), parameter=t)
length_r_gamma = gamma.length(interval=interval, manifold=surfaceR)
print(length_r_gamma)

# Using formula 1.55
"""
L_E(r(phi_inv(eta))) = int(sqrt(eta'(t), G_v(eta(t)), eta'(t))))
"""
y1, y2, y3= sp.symbols('y1, y2, y3')
ySpace = Space(metric=sp.Matrix([[4*y1**2 + 1, 0], [0, y1**2]]), 
               coord_vars = sp.Matrix([y1, y2]))

eta = Curve(curve_expr=sp.Matrix([1, t]), parameter=t)
length_r_phiinv_eta = eta.length(interval=interval, manifold=ySpace)
print(length_r_phiinv_eta)

# both are equal to 2*pi
