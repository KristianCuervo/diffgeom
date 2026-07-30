import diffgeom as dg
import sympy as sp

# Exam 2019 Problem 4

x1, x2 = sp.symbols('x1, x2')
space = dg.Space(metric=(x1**2)*sp.eye(2), 
                 coord_vars=[x1, x2])

# 4.1 find coordinate function R1221(x1, x2) of 
# curvature tensor R for the given metric tensor field g

R1221 = space.riemanns[0, 1, 1, 0]
print(R1221)

# 4.2 find the sectional curavture k 
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
K = space.sectional_curvature(X=e1, Y=e2)
print(K)
