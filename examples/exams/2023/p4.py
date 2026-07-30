import diffgeom as dg
import sympy as sp

# Exam 2023 Problem 4

t = sp.symbols('t')
gamma1 = sp.Matrix([t, 0])
gamma2 = sp.Matrix([0, t])

chris111 = 1
chris122 = 1

"""
Looking at def 3.22 and eq 3.49, we find that the 
second derivative of gamma1 and gamma2 is zero, and that the sum 
is only non-zero for i,j=1 and i,j=2. We then use this to
do the sum of gammaprime[i] * gammaprime[j] for  these values
"""
acc_gamma1 = sp.Matrix([gamma1.diff(t)[0] * gamma1.diff(t)[0] + gamma1.diff(t)[1]*gamma1.diff(t)[1], 0 ])
acc_gamma2 = sp.Matrix([gamma2.diff(t)[0] * gamma2.diff(t)[0] +  gamma2.diff(t)[1]*gamma2.diff(t)[1], 0 ])
# look at k=1 for both, thus the sums are both in k=1
print(acc_gamma1, acc_gamma2) # [1, 0], [1, 0]

# 3.2 Which can be reparamtereised as a geodesic
"""
We try to see whether the curves as pre-geodesics, by looking at whether
we can satisfy a rho(t) which is a smooth function such that:
D/dt(gammaprime) = rho(t) * gammaprime(t)
"""
print(acc_gamma1 == 1*gamma1.diff(t))
"""
We can actually observe that for rho(t)=1, that gamma1 is a pre-geodesic 
pre-geodesics. Thus they can be reparamterised into geodesic curves.
However for 

[1, 0] = rho(t) * [0, 1]
which is the case for gamma2, then there is clearly no case where this can
be reparameterised to another solution.
"""
