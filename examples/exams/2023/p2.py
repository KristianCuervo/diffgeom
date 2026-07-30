import diffgeom as dg
import sympy as sp

# 2023 Problem 2

x1, x2 = sp.symbols('x1 x2') 
y1, y2 = sp.symbols('y1 y2')

euclid2d = dg.Space(metric=sp.eye(2),
                    coord_vars=sp.Matrix([x1, x2]))

dummy_space = dg.Space(metric=sp.eye(2),
                       coord_vars=sp.Matrix([y1, y2]))

phi = dg.Map(expr=sp.Matrix([sp.sin(x1), x2]),
             domain_space=euclid2d, 
             codomain_space=dummy_space)

space_V = phi.pushforward_metric()
print(space_V.metric.subs({y1: x1,
                           y2: x2})) # Replace from dummy variables back to first variables
