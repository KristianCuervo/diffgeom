import diffgeom as dg
import sympy as sp

# Exam 2020 Problem 1

x1, x2 = sp.symbols('x1 x2')
k = sp.symbols('k', real=True, positive=True)

G = sp.Matrix([[1, 0],
               [0, k**2 + x1**2]])

space = dg.Space(metric=G,
                 coord_vars=sp.Matrix([x1, x2]))
p = (0, 0)
I, It, t = space.find_indicatrix(point=p)
print(I)
parameterisation = sp.Matrix([sp.cos(t), (1/k)*sp.sin(t)]) # parameterisation'
print("The parameterisation is equal to: ", parameterisation)
test = space.metric_tensor(V=parameterisation, W=parameterisation, metric=space.metric.subs({x1:0,
                                                                                             x2:0}))
print(sp.trigsimp(test[0])) # parameterisation holds and is equal to one


# 1.2 Show that e1 is not a killing field
e1 = dg.Field(vector_field=sp.Matrix([1, 0]),
              space=space)
print(e1.is_killing())

# 1.3 Show that e2 IS a killing field
e2 = dg.Field(vector_field=sp.Matrix([0, 1]),
              space=space)
print(e2.is_killing())

# 1.4 let f=x1+ x2. Find gradient of f
f = dg.Field(vector_field=sp.Matrix([x1 + x2]),
             space=space)
print(f.gradient)

# 1.5 Show non-zero christoffel symbols
print(space.christoffels[0,1], space.christoffels[1, 0])
print(space.christoffels[1,1])


# 1.6 Show that curve is a geodeseic
s = sp.symbols('s', real=True)
alpha = sp.symbols('alpha', real=True)
gamma = dg.Curve(curve_expr=sp.Matrix([s, alpha]),
                 parameter=s,
                 manifold=space)
print(gamma.is_geodesic()) # truee

# new curve
t = sp.symbols('t', real=True)
beta = sp.symbols('beta', real=True)
mu = dg.Curve(curve_expr=sp.Matrix([beta, t]),
              parameter=t,
              manifold=space)
print(mu.speed)
# parameterisation: s = t * mu.speed, t = s/mu.speed
eta = dg.Curve(curve_expr=sp.Matrix([beta, s/mu.speed]),
                                    parameter=s,
                                    manifold=space)
print(eta.speed) # equals to 1

# 1.8 Find acceleration of eta

acc_eta = eta.acceleration_vector()
print(acc_eta)

# 1.9 is a riemannian circle
kappa = space.metric_tensor(V=acc_eta,
                            W=acc_eta,
                            metric=space.metric.subs({
                                space.vars[0]: eta.expr[0],
                                space.vars[1]: eta.expr[1]
                            }))

print(kappa)
# which is a constant non-zero beta for all beta as it is not dependent on variables
# 1.10 is a geodesic for one value of beta
"""
Only when acceleration is equal to zero, which is when beta=0
"""
print(kappa.subs(beta, 0)) # 0

# 1.11 Find coordinate function r1221

R1221 = space.riemanns[0,1,1,0]
print(R1221)

# 1.12 Find sectional curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
K = space.sectional_curvature(X=e1, Y=e2)
print(K)

# Show that K--> 0 as x1 --> +- inf
