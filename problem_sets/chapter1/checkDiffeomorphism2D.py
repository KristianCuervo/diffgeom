import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

"""
Exercise 1.13 
"""

X = sp.Matrix([sp.symbols('x1'), sp.symbols('x2')])
phi = dg.Map(sp.Matrix([sp.exp(X[0]) + X[1], X[0]]), X)

# a) show that phi is a diffeomorphism
print(phi.isDiffeomorphism)

# b) Find Jphi and Jphi^-1
print("Jphi = ", phi.j)

finv, Y = phi.inverseMapping()
phiInv = dg.Map(finv, Y)
print("Jphi^-1 = ", phiInv.j)

# c) Show that these jacobian matrices are the inverses of each other
print("J^-1 phi = ", phi.j.inv())
print("JphiInv^in X", phiInv.j.subs({Y[0]: phi.expr[0], Y[1]: phi.expr[1]}))

if phi.j.inv() == phiInv.j.subs({Y[0]: phi.expr[0], Y[1]: phi.expr[1]}):
    print("The jacobians are inverses of each other")