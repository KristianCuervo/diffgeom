import diffgeom as dg
import sympy as sp

"""
Exercise 1.13

phi(x1, x2) = (e^x1 + x2, x1) maps R2 to R2. Show it is a diffeomorphism,
find its Jacobian and the Jacobian of its inverse, and confirm that the two
are inverses of each other.
"""

x1, x2 = sp.symbols('x1 x2')
X = sp.Matrix([x1, x2])

y1, y2 = sp.symbols('y1 y2')
Y = sp.Matrix([y1, y2])

# A Map needs a domain and a codomain. Neither metric matters here -- only
# the Jacobians do -- so both are flat R2.
domain = dg.Space(metric=sp.eye(2), coord_vars=X)
codomain = dg.Space(metric=sp.eye(2), coord_vars=Y)

phi = dg.Map(expr=sp.Matrix([sp.exp(x1) + x2, x1]),
             domain_space=domain,
             codomain_space=codomain)

# a) show that phi is a diffeomorphism
print("phi is a diffeomorphism:", phi.is_diffeomorphism)

# b) Find Jphi and Jphi^-1
print("Jphi =", phi.jacobian_matrix)

phi_inv = phi.get_inverse_map()
print("phi^-1 =", phi_inv.expr)
print("Jphi^-1 =", phi_inv.jacobian_matrix)

# c) Show that these jacobian matrices are the inverses of each other.
# Jphi^-1 is expressed in the codomain coordinates y, so pull it back to x
# by substituting y = phi(x) before comparing.
back_in_x = phi_inv.jacobian_matrix.subs({y1: phi.expr[0], y2: phi.expr[1]})
print("(Jphi)^-1  =", sp.simplify(phi.jacobian_matrix.inv()))
print("Jphi^-1|_x =", sp.simplify(back_in_x))

if sp.simplify(phi.jacobian_matrix.inv() - back_in_x) == sp.zeros(2, 2):
    print("The jacobians are inverses of each other")
