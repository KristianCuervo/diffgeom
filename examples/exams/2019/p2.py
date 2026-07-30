import diffgeom as dg
import sympy as sp

# Exam 2019 Problem 2

s = sp.symbols('s')
L = sp.symbols('L')
arc_length = (0, L)


fs = s + s**3
jacobi = dg.JacobiField(f_s=fs, parameter=s)
ks = jacobi.find_sectional_curvature().rhs
ks_0 = ks.subs(s, 0)
"""
The sectional curvature at point s=0 is given by:
"""
print(ks_0)
