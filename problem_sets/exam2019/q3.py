import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2019 Problem 3
"""
Poincare half plane model
"""
x1 = sp.symbols('x1')
x2 = sp.symbols('x2', positive=True)

space = dg.Space(metric=(1/x2**2)*sp.eye(2),
                 coord_vars=sp.Matrix([x1, x2]))


# 3.1 Show that V is a killing field in M2
V = dg.Field(vector_field=sp.Matrix([1, 0]),
             space=space)
V.is_killing()

# 3.2 Show that W is not a killing field
W = dg.Field(vector_field=sp.Matrix([0, 1]),
             space=space)
W.is_killing()

# 3.3 Reparameterise curve C
"""
The curve is given by x1**2 + x2**2 = 1
This is the equation of a circle (but it half-circle as x2>0)
We can reparameterise this to a curve as 
x1 = cos(t), x2=sin(t), t=[0, pi]
"""
t = sp.symbols('t')
parameterisation = sp.Matrix([sp.cos(t), sp.sin(t)])
gamma = dg.Curve(curve_expr=sp.Matrix([sp.cos(t), sp.sin(t)]),
                 parameter=t,
                 manifold=space)

"""
A(x1, x2) = [-x2, x1] --> A(t) = [-sin(t), cos(t)]
B(x1, x2) = [x1, x2] -->  B(t) = [cos(t), sin(t)]
"""
A_t = sp.Matrix([-sp.sin(t), sp.cos(t)])
B_t = sp.Matrix([sp.cos(t), sp.sin(t)])

f = sp.Function('f')(t)
h = sp.Function('h')(t)

covar_A = space.covariant_derivative(vector=f*A_t, curve=gamma)
eqs = [sp.Eq(covar_A[0], 0),
       sp.Eq(covar_A[1], 0)]
print(eqs[0])
sol1 = sp.dsolve(eqs[0], f)
sol2 = sp.dsolve(eqs[1], f)
print(sol1, sol2)
print(sol1 == sol2)
# f is equation to c1*sin(t)

# 3.5
covar_B = space.covariant_derivative(vector=h*B_t, curve=gamma)
eqs = [sp.Eq(covar_B[0], 0),
       sp.Eq(covar_B[1], 0)]
sol3 = sp.dsolve(eqs[0], h)
sol4 = sp.dsolve(eqs[1], h)
print(sol3, sol4)
print(sol3 == sol4)
# h is equation to c2*cos(t)


