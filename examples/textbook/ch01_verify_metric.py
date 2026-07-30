import diffgeom as dg
import sympy as sp

"""
Exercise 1.23

The paraboloid r(x1, x2) = (x1, x2, x1^2 + x2^2) sits in Euclidean R3.
Check that the metric it inherits from the ambient space matches the
metric given in the exercise.
"""

x1, x2, x3 = sp.symbols('x1 x2 x3')
X = sp.Matrix([x1, x2])

G_given = sp.Matrix([[4*x1**2 + 1, 4*x1*x2],
                     [4*x1*x2, 4*x2**2 + 1]])
print("G_given:")
print(G_given)

euclid3d = dg.Space(metric=sp.eye(3), coord_vars=sp.Matrix([x1, x2, x3]))
S = dg.Surface(embedded_expr=sp.Matrix([x1, x2, x1**2 + x2**2]),
               surface_vars=X,
               ambient_space=euclid3d)
G_computed = S.metric
print("G_computed:")
print(G_computed)

if G_given.equals(G_computed):
    print("The metric tensors are equal.")
