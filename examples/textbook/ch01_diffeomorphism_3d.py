import diffgeom as dg
import sympy as sp

"""
Exercise 1.15

phi(x1, x2, x3) = (x2, x3, x1 + x2 + x3) is a linear map on R3. Show it is a
diffeomorphism, find its Jacobian and the Jacobian of its inverse, and confirm
that the two are inverses of each other.
"""

x1, x2, x3 = sp.symbols('x1 x2 x3')
X = sp.Matrix([x1, x2, x3])

y1, y2, y3 = sp.symbols('y1 y2 y3')
Y = sp.Matrix([y1, y2, y3])

domain = dg.Space(metric=sp.eye(3), coord_vars=X)
codomain = dg.Space(metric=sp.eye(3), coord_vars=Y)

phi = dg.Map(expr=sp.Matrix([x2, x3, x1 + x2 + x3]),
             domain_space=domain,
             codomain_space=codomain)

# a) show that phi is a diffeomorphism
print("phi is a diffeomorphism:", phi.is_diffeomorphism)

# b) find Jphi and Jphi^-1
print("Jphi =", phi.jacobian_matrix)

phi_inv = phi.get_inverse_map()
print("phi^-1 =", phi_inv.expr)
print("Jphi^-1 =", phi_inv.jacobian_matrix)

# c) Show that these jacobian matrices are the inverses of each other.
# phi is linear, so its Jacobian is constant and no pullback is needed --
# but substituting y = phi(x) is harmless and matches the 2D case.
back_in_x = phi_inv.jacobian_matrix.subs({y1: phi.expr[0],
                                          y2: phi.expr[1],
                                          y3: phi.expr[2]})
print("(Jphi)^-1  =", sp.simplify(phi.jacobian_matrix.inv()))
print("Jphi^-1|_x =", sp.simplify(back_in_x))

if sp.simplify(phi.jacobian_matrix.inv() - back_in_x) == sp.zeros(3, 3):
    print("The jacobians are inverses of each other")
