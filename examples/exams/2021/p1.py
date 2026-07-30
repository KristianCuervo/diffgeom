import diffgeom as dg
import sympy as sp

# Exam 2021 Problem 1

x1, x2 = sp.symbols('x1, x2')
G = sp.Matrix([[1, 0],
               [0, (1+x1**2)**2]])
space = dg.Space(metric=G, 
                 coord_vars=sp.Matrix([x1, x2]))

# 1.1 Find parameterisation of I(1,0) for g s.t
# gamma(t) = f(t)e1 + h(t)e2 tin[-pi,pi]
p = (1,0)
I, It, t = space.find_indicatrix(point=p)
parameterisation = sp.Matrix([sp.cos(t), sp.sin(t)])
test = space.metric_tensor(V=parameterisation,
                           W=parameterisation,
                           metric=space.metric.subs({x1:p[1],
                                                    x2:p[1]}))
print(sp.trigsimp(test)) # 1 --> Thus a parameterisation which holds for the indicatrix


# 1.2 Find all non-zero christoffel symbols
print(space.christoffels[0,1], space.christoffels[1,0])
print(space.christoffels[1,1])
"""
chris12, = chris21 = 2.0*x1/(x1**2 + 1)
chris22 = -2.0*x1*(x1**2 + 1)
"""

# 1.3 Show that eta is a unit speed parameterised curve and not a geodesic

s = sp.symbols('s', real=True)
eta = dg.Curve(curve_expr=sp.Matrix([1, s/2]),
               parameter=s,
               manifold=space)
print(eta.speed) # 1
print(eta.is_geodesic()) # False

# 1.4 Find all possible coordinate functions s.t V is parallel along eta
"""
V(s) = v1(s), v2(s), is a vector field along eta
For it to be a parallel function, its covariant derivative must 
be equal to zero.
"""
v1 = sp.Function('v1')(s)
v2 = sp.Function('v2')(s)
V = sp.Matrix([v1, v2])
covar_V = dg.Field(vector_field=V, space=space).covariant_derivative(curve=eta)
#print(covar_V)
eqs = [sp.Eq(covar_V[i], 0) for i in range(space.dim)]
sol = sp.dsolve(eqs, [v1, v2])
v1_sol = sol[0].rhs
v2_sol = sol[1].rhs

"""
Now we verify that this is a parallel transport field by yet
again making sure its covariant derivative is zero
"""
V_sol = sp.Matrix([v1_sol, v2_sol])
parallel_test = dg.Field(vector_field=V_sol, space=space).covariant_derivative(curve=eta)
print(parallel_test) # 0 ! 


# 1.5 Find all pairs of functions f and h s.t Z is a killing field
"""
Z(x1,x2)=(f(x2), h(x1))
"""
f = sp.Function('f')(x2)
h = sp.Function('h')(x1)
Z = sp.Matrix([f, h])
Zfield = dg.Field(vector_field=Z,
                  space=space)
Zkilling = Zfield.is_killing(result=True)
for condition in Zkilling:
    print(condition," = 0")
"""
Read-off:
For this to be true. It must be such that firstly, f(x2)=0.
Then given this we solve for h1_derivative =  (x1**2 + 1)**(-2)
"""
hprime = h.diff(x1)*(x1**2 + 1)**(-2)
h_sol = sp.dsolve(sp.Eq(hprime,0), h)
print(h_sol) # h = c1, any constant. Then hprime is 0.

# 1.6 find K
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])
K = space.sectional_curvature(X=e1, Y=e2)
print(K)
"""
The order of the numerator is x1**2, and the order 
of the denominator is x1**4. As x1-->+-inf, then you have
that (1/x**2)--> 0. 

"""
