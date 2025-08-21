import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Problem 9.18

y1, y2, y3 = sp.symbols('y1, y2, y3')
x1, x2 = sp.symbols('x1, x2')
R = sp.symbols('R', real=True, positive=True)

ambient_space = dg.Space(metric=(1/(y3**2))*sp.eye(3), 
                                 coord_vars=[y1, y2, y3])
plane_space = dg.Space(metric=sp.eye(2),
                       coord_vars=[x1, x2])

theta =sp.symbols('theta')
phi = dg.Map(expr=sp.Matrix([x1, 
                             x2*sp.cos(theta),
                             1+x2*sp.sin(theta)]),
            domain_space=plane_space,
            codomain_space=ambient_space)

Gu = plane_space.metric_tensor(V=phi.jacobian_matrix, W=phi.jacobian_matrix,
                               metric=ambient_space.metric.subs({y1: phi.expr[0],
                                                                y2: phi.expr[1],
                                                                y3: phi.expr[2]}))
print(sp.trigsimp(Gu)) # BAM


# 9.22

phi2 = dg.Map(expr=sp.Matrix([R*sp.cos(x2),
                              R*sp.sin(x2),
                              x1]),
              domain_space=plane_space,
              codomain_space=ambient_space)
# Calculate the intrinsic metric
Gu2 = plane_space.metric_tensor(V=phi2.jacobian_matrix, W=phi2.jacobian_matrix,
                               metric=ambient_space.metric.subs({y1: phi2.expr[0],
                                                                y2: phi2.expr[1],
                                                                y3: phi2.expr[2]}))
print(sp.trigsimp(Gu2))


# calculate instrinsic curvature
e1 = sp.Matrix([1, 0])
e2 = sp.Matrix([0, 1])

induced_space = dg.Space(metric=Gu2, 
                         coord_vars=[x1, x2])

sec = induced_space.sectional_curvature(X=e1, Y=e2)
print(sp.trigsimp(sec))



# try with new surface

surface_phi = dg.Surface(embedded_expr=sp.Matrix([[R*sp.cos(x2),
                              R*sp.sin(x2),
                              x1]]),
                        surface_vars=[x1, x2],
                        ambient_space=ambient_space)

print("GU with surface class", surface_phi.metric) # same GU as before
print("Sectional Curvature w/ surface class", surface_phi.sectional_curvature(X=e1, Y=e2)) # same as before

