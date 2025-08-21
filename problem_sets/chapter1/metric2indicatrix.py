import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exercise 1.53

t = sp.symbols('t')
u = sp.Matrix([sp.cos(t), sp.sin(t)])
T1 = sp.Matrix([[3, -7*sp.sqrt(3)],
               [3*sp.sqrt(3), 7]])

G1 = sp.simplify(T1.T.inv() * T1.inv())
print(G1)

T_homework = sp.Matrix([[-sp.Rational(1, 4), sp.sqrt(3)/2],
                        [sp.sqrt(3)/2, sp.Rational(1, 2)]])
print(T_homework.inv())
G2 = sp.simplify(T_homework.T.inv() * T_homework.inv())
print(G2)