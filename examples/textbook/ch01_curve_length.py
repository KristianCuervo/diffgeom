from diffgeom import Surface, Curve, Space, Map
import sympy as sp

"""
Exercise 1.29

The paraboloid r(x1, x2) = (x1, x2, x1^2 + x2^2) in Euclidean R3.
Measure the same curve two ways and check the answers agree:

  formula 1.48 -- integrate the g-speed in the surface chart (x1, x2)
  formula 1.55 -- integrate in polar coordinates (y1, y2) instead

Both give 2*pi.
"""

x1, x2, x3 = sp.symbols('x1 x2 x3')
euclid3d = Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))

# The surface, and the metric it inherits from R3
surfaceR = Surface(embedded_expr=sp.Matrix([x1, x2, x1**2 + x2**2]),
                   surface_vars=sp.Matrix([x1, x2]),
                   ambient_space=euclid3d)
print(surfaceR.metric)  # a) Gu is correct

# Polar coordinates on the chart: phi(x1, x2) = (r, theta)
y1, y2 = sp.symbols('y1 y2')
ySpace = Space(metric=sp.Matrix([[4*y1**2 + 1, 0],
                                 [0, y1**2]]),
               coord_vars=sp.Matrix([y1, y2]))

phiMap = Map(expr=sp.Matrix([sp.sqrt(x1**2 + x2**2), sp.arg(x1 + x2*sp.I)]),
             domain_space=surfaceR,
             codomain_space=ySpace)
print("phi is a diffeomorphism:", phiMap.is_diffeomorphism)

interval = (-sp.pi, +sp.pi)
t = sp.symbols('t')

# Using formula 1.48: the unit circle in the (x1, x2) chart
gamma = Curve(curve_expr=sp.Matrix([sp.cos(t), sp.sin(t)]),
              parameter=t,
              manifold=surfaceR)
length_r_gamma = gamma.length(interval=interval, manifold=surfaceR)
print(length_r_gamma)

# Using formula 1.55
"""
L_E(r(phi_inv(eta))) = int(sqrt(eta'(t), G_v(eta(t)), eta'(t))))
"""
eta = Curve(curve_expr=sp.Matrix([1, t]),
            parameter=t,
            manifold=ySpace)
length_r_phiinv_eta = eta.length(interval=interval, manifold=ySpace)
print(length_r_phiinv_eta)

# both are equal to 2*pi
