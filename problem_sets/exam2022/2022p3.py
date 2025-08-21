import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2022 Problem 3
"""
Solved by linearity. Look if T(fX, Y, Z) == f T(X, Y, Z). which it is not.
"""

