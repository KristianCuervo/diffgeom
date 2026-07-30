import diffgeom as dg
import sympy as sp

import numpy as np
# Exam 2022 Problem 2
# Which matrices COULD be the metric matrix function
"""
The metric matrix function must be:
symmetric
positive definite.

G1 is clearly not symmetric. G2 and G3 are symmetric but 
we look for if they are positive definie
"""

# G2:
"""
G2 can take negative values in the diagonal. Thus it can 
form negative eigenvalues based on. We test this with (x1, x2) = (-101, 11)
"""
test_g2 = np.array([[-10, 1],
                    [1, -10]])
print(np.linalg.eig(test_g2)[0]) # this gives negative eigenvalues
"""
Thus G2 does not form an LRM.
"""

# G3: 
"""
G3 does form an LRM. This is because it is symmetric, and the 
off-diagonal values are limited by the cosine function, which
has a minimum of (0) and maximum of (1). Given that all
values in the matrix will take positive eigenvalues at all times,
then the matrix will only yield positive eigenvalues, and 
is positive definite.
"""
