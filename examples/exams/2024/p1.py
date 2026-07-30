import diffgeom as dg
import sympy as sp

# Exam 2024 Problem 1

x1, x2, x3 = sp.symbols('x1 x2 x3')
euclid3d = dg.Space(metric=sp.eye(3),
                    coord_vars=sp.Matrix([x1,x2,x3]))
surface_r = dg.Surface(embedded_expr=sp.Matrix([sp.cos(x2)*sp.cos(x1), sp.cos(x2)*sp.sin(x1), sp.sin(x2)]),
                       surface_vars=sp.Matrix([x1, x2]),
                       ambient_space=euclid3d)
# 1.1 Metric of surface
print(surface_r.metric) # Holds

# 1.2
"""
The metric is only dependent on x2. Thus in the direction 
of e1=[1,0], then a vector field would be a killing field 
as its vectors are preserved in length along the direction of path.

So the given three vector fields, the only one
which points in the e1 direction is [7, 0]. This would be a killing field.
"""
test_field = dg.Field(vector_field=sp.Matrix([7, 0]),
                      space=surface_r)
print(test_field.is_killing()) # True

# 1.3 Non-zero christoffel symbols
print(surface_r.christoffels)

# 1.4 Shortest curve
t = sp.symbols('t')
a = sp.symbols('a')
gamma = dg.Curve(curve_expr=sp.Matrix([t, a]),
                 parameter=t,
                 manifold=surface_r)

acc_gamma = gamma.acceleration_vector()
print(acc_gamma)
"""
If acc(gamma) == 0, then gamma is a geodesic. Thus we look at for
what values of a this is true.

As acc = [[0], [0.5*sin(2*a)]]. This holds for a = 0.
"""
# 3.5 Find coordinate function R1221 in U
R1221 = surface_r.riemanns[0,1,1,0]
print(R1221) # 1

# 3.6 
field_W = dg.Field(vector_field=sp.Matrix([sp.sin(t*sp.sin(a))/sp.cos(a),
                                           sp.cos(t*sp.sin(a))]),       
                    space=surface_r)

covar_W = field_W.covariant_derivative(curve=gamma)
print(covar_W) # 0
"""
The covariant derivative is zero and thus the choice
of alpha does not affect whether the field W is a 
parallel transport field or not, it will always be a
parallel transport field along gamma. 
"""

# 3.7 
y1, y2 = sp.symbols('y1 y2')
dummy_space = dg.Space(metric=sp.Matrix([[sp.cos(x2)**2, 0],
                                         [0, 1]]),
                       coord_vars=sp.Matrix([y1, y2]))
space_U = dg.Space(metric=sp.Matrix([[sp.cos(x2)**2, 0],
                                         [0, 1]]),
                   coord_vars=sp.Matrix([x1, x2]))
phi = dg.Map(expr=sp.Matrix([x1, sp.ln(sp.tan(sp.pi/4 + x2/2))]),
             domain_space=space_U,
             codomain_space=dummy_space)
print(phi.expr)
phi_inv = phi.get_inverse_map()
print(phi_inv.expr)

another_V = phi.fromUtoV()
print(sp.trigsimp(another_V.metric))

# 3.9
"""
{7,0} as it is independent of y1}
"""

# 3.10 
R1221 = another_V.riemanns[0, 1, 1, 0]
print(R1221)
