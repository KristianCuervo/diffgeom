import diffgeom as dg
import sympy as sp

# Exercise 9.23

y1, y2, y3 = sp.symbols('y1, y2, y3')
ambient_space = dg.Space(sp.Matrix([[1+y1**2, 0, 0],
                                    [0, 1, 0],
                                    [0, 0, 1]]),
                        coord_vars=[y1, y2, y3])

x1, x2 = sp.symbols('x1, x2')
plane_space = dg.Space(metric=sp.eye(2),
                       coord_vars=[x1, x2])
a, b = sp.symbols('a, b')
phi = dg.Map(sp.Matrix([x1, x2, sp.Rational(1, 2)*(a*x1**2 + b*x2**2)]),
             domain_space=plane_space,
             codomain_space=ambient_space)

Gu = plane_space.metric_tensor(V=phi.jacobian_matrix, W=phi.jacobian_matrix,
                               metric=ambient_space.metric.subs({y1: phi.expr[0],
                                                                 y2: phi.expr[1],
                                                                 y3: phi.expr[2]}))
print(Gu) # Correct

# find the intrinsic curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])

induced_space = dg.Space(metric=Gu, coord_vars=[x1, x2])
sec = induced_space.sectional_curvature(X=e1, Y=e2)
print(sp.simplify(sec))  # correct, tricky with which simplify to use
