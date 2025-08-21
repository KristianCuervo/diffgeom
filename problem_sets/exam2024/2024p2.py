import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2024 Problem 1

s = sp.symbols('s')
fs = s*(1+sp.sin(s)**2)**4
jacob = dg.JacobiField(f_s=fs, parameter=s)
print(sp.simplify(sp.trigsimp(jacob.find_sectional_curvature())))



