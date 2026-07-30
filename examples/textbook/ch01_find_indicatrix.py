import diffgeom as dg
import sympy as sp

"""
Find the indicatrix of a constant metric.

The indicatrix at a point p is the set of vectors V with I(p)(V) = 1, where
I(p)(V) = V^T G(p) V. find_indicatrix returns that quadratic form, the same
form restricted to the unit circle V = (cos t, sin t), and the parameter t.

Both metrics below are constant, so the indicatrix is the same at every
point and the choice of p is arbitrary.
"""

x1, x2 = sp.symbols('x1 x2')
coords = sp.Matrix([x1, x2])
origin = (0, 0)

Gu = sp.Matrix([[6, 2],
                [2, 1]])
space = dg.Space(metric=Gu, coord_vars=coords)
Ip, Ipt, t = space.find_indicatrix(point=origin)
print("I(p)(V)     =", sp.simplify(Ip))
print("on V=(cos t, sin t):", sp.simplify(Ipt))

homework_G = sp.Matrix([[8, sp.sqrt(3)],
                        [sp.sqrt(3), 6]])
homework_space = dg.Space(metric=homework_G, coord_vars=coords)
Ip2, Ipt2, t2 = homework_space.find_indicatrix(point=origin)
print("I(p)(V)     =", sp.simplify(Ip2))
print("on V=(cos t, sin t):", sp.simplify(Ipt2))
