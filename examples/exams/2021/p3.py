import diffgeom as dg
import sympy as sp

# Exam 2021 Problem 3
x1 = sp.symbols('x1')
x2 = sp.symbols('x2', positive=True)
euclid = dg.Space(metric=sp.eye(2),
                  coord_vars=sp.Matrix([x1,x2]))

y1, y2 = sp.symbols('y1, y2')
image_euclid = dg.Space(metric=sp.eye(2),
                        coord_vars=sp.Matrix([y1, y2]))
# 3.1 Describe the image 
"""
Given that x1 can be any real number, and x2 is strictly positive. 
Then if y1 = 2x1 + x2, and y2 = x1 + x2, we can actually allow
y1 and y2 to be anywhere in R.
"""
# 3.2 Check diffeo and find inverse expression
phi = dg.Map(expr=sp.Matrix([2*x1 + x2, x1 + x2]),
             domain_space=euclid,
             codomain_space=image_euclid)
print(phi.is_diffeomorphism) # True

phi_y = phi.get_inverse_map()
print(phi_y.expr) # Inverse expression 

# 3.3 Find jacobian in terms of phi_inv
print(phi_y.jacobian_matrix)

# 3.4 What is the metric in Gv(y1, y2)
Gu = ((1/x2)**2)*sp.eye(2)
space_u = dg.Space(metric=Gu,
                   coord_vars=sp.Matrix([x1, x2]))
phi2 = dg.Map(expr=sp.Matrix([2*x1 + x2, x1 + x2]),
              domain_space=space_u,
              codomain_space=image_euclid)
space_v = phi2.pushforward_metric()
print(space_v.metric)

# 3.5 Find the parameterisation of the geodesics of U, in V.

s = sp.symbols('s', real=True)
B = sp.symbols('B', real=True)
eta = sp.Matrix([B, sp.exp(s)])
gamma = phi2.expr.subs({x1: eta[0],
                        x2: eta[1]})
print(gamma) # B = 1
print(gamma.diff(s)) # Meets gamma'(0) = 1,1
gamma_sol = sp.Matrix([2 + sp.exp(2), 1 + sp.exp(s)]) # Answer 3.5

# 3.6 Same but for mu
C = sp.symbols('C')
K = sp.symbols('K')
xi = sp.Matrix([C*sp.tanh(s) + K, C/sp.cosh(s)])
mu = phi2.expr.subs({x1: xi[0],
                     x2: xi[1]})
print(mu.subs(s, 0)) # [C+2K, C+K]
"""
To solve the above equation, know that C+2k=3, and C+k=2
Then C=1, K=1 solves this.
Making sure that the derivative condition is also met
"""
print(mu.diff(s).subs(s,0)) # [2c, c] 
"""
This condition is also met for c=1
"""
