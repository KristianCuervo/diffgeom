import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# 2022 Problem 1 #
x1,  x2 = sp.symbols('x1 x2')
X = sp.Matrix([x1, x2])
k = sp.symbols('k', real=True, positive=True)

Gk = sp.Matrix([[1, 0], [0, k**2 + x1**2]])

space = dg.Space(metric=Gk, coord_vars=X)

# ----- 1.4 ------ #
# Find gradient of f

f = dg.Field(vector_field=sp.Matrix([x1 + x2]), space=space)
print(f.gradient)

# ----- 1.5 ------ #
# Show non-zero christoffel symbols
print(space.christoffels)

# ----- 1.6 ------ #
# Show curve (gamma(s) = (s, alpha) is a geodesic for all alpha)
alpha = sp.symbols('alpha', real=True)
s = sp.symbols('s', real=True)
gamma = dg.Curve(curve_expr = sp.Matrix([s, alpha]), parameter=s)

print(gamma.is_geodesic(space)) 

# ----- 1.7 ------ #
# Reparameters mu to become arc length parameterised (eta)
beta = sp.symbols('beta', real=True)
t = sp.symbols('t', real=True)
mu_beta = dg.Curve(curve_expr=sp.Matrix([beta, t]), parameter=t)

mu_speed = sp.sqrt(space.metric_tensor(V=mu_beta.derivative,
                               W=mu_beta.derivative,
                               metric=space.metric.subs({
                                    space.vars[0]: mu_beta.expr[0],
                                    space.vars[1]: mu_beta.expr[1]
                               })))
print(mu_speed) #beta**2 + k**2
# s = (beta**2 + k**2) * t
# then t = s / (beta**2 + k**2)
s = sp.symbols('s', real=True)
eta = dg.Curve(curve_expr=sp.Matrix([beta, s/sp.sqrt((beta**2 + k**2))]), parameter=s)
eta_speed = space.metric_tensor(V=eta.derivative,
                                W=eta.derivative,
                                metric=space.metric.subs({
                                    space.vars[0]: eta.expr[0],
                                    space.vars[1]: eta.expr[1]
                                }))
print(eta_speed) # 1 --> now arc length parameterised

# ----- 1.8 ------ #
# Find acceleration of eta for all s and all beta
acc_eta = space.acceleration_vector(curve=eta)
print(acc_eta) 

# ----- 1.9 ------ #
# Show eta is a riemannian circle for beta =! 0
"""
From def 5.8: an arc length parameterised riemannian circle in an LRM
is an arc length parameterised Riemmanian helix with torsion tau=0 
and constant curvature kappa.

We look at the curvature kappa: which is the norm of the acceleration.
"""
kappa = space.metric_tensor(V=acc_eta,
                            W=acc_eta,
                            metric=space.metric.subs({
                                space.vars[0]: eta.expr[0],
                                space.vars[1]: eta.expr[1]
                            }))
print(kappa)
"""
kappa = 1.0*beta**2/(beta**2 + k**2)**2
which for non-zero beta is a constant curvature.
"""

# ----- 1.10 ----- #
# Show that eta is a geodesic for one value of beta
"""
eta is a geodesic is the acceleration is zero.
given that the accerlation is:
Matrix([[-1.0*beta/(beta**2 + k**2)], [0]])

this is only equal to zero if beta = 0.
"""

# ----- 1.11 ----- #
# Find coordinate function R1221(x1, x2) of curavture tensor R for given metric gk
R1221 = space.riemanns[0,1,1,0]
print(R1221) #-1.0*k**2/(k**2 + x1**2)

# ----- 1.12 ----- #'
#Find the sectional curvature function K for gk
K = space.sectional_curvature(X=sp.Matrix([1, 0]), Y=sp.Matrix([0, 1]), R=R1221)
print(K) #1.0*k**2/(k**2 + x1**2)**2

# ----- 1.13 ----- #
# Show limit goes to 0 for x1--> +-inf
"""
Given that the sectional curvature is: 1.0*k**2/(k**2 + x1**2)**2
the lim(x1->+-inf) K = 0 as the denominator 
will continuously grow and shrink the fraction to 0.
"""
