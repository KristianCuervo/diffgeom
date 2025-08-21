import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp


x1, x2 = sp.symbols('x1, x2')
Gu = sp.Matrix([[6, 2], [2, 1]])
space = dg.Space(metric=Gu, coord_vars=sp.Matrix([x1, x2]))
i, t = space.find_indicatrix()
print(sp.simplify(i))


homework_G = sp.Matrix([[8, sp.sqrt(3)], [sp.sqrt(3), 6]])
newThing = dg.Space(metric=homework_G, coord_vars=sp.Matrix([x1, x2]))
i2, t2 = newThing.find_indicatrix()
print(sp.simplify(i2))
