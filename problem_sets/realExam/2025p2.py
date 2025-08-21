import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import differentialgeometry as dg
import sympy as sp

# Exam 2025 Problem 2


x1 = sp.symbols('x1', real=True, positive=True)
x2 = sp.symbols('x2', real=True)

space = dg.Space(metric=sp.Matrix([[2, 0],
                                   [0, x1**2]]),
                coord_vars=sp.Matrix([x1, x2]))

# 3.1 Find all vector fields which are a killing field in this space

f_x1 = sp.Function('f')(x1)
h_x2 = sp.Function('h')(x2)
w_field = dg.Field(vector_field=sp.Matrix([f_x1, h_x2]),
                   space=space)
killing_test = w_field.is_killing(result=True)
print(sum(killing_test))
print(killing_test)
"""
We require that the above sum(killing_test) is equal to zero. This is followed
from the construction given in proposition 2.38. 

From observation we can state the following:
1) 8 * f'(x1) = 0 --> This requires f'(x1)=0

2) 4*x1*f(x1) + 2*x1**2 * h'(x2) = 0
# from 1) we know that f(x1) = constant. Thus we can't match
the terms x1 and x1**2 in the two terms. Thus we also require that
f(x1) = 0 and h'(x2) = 0 

This leads us with
w_field = [0, c] for any constant c in Real
We test this with the following:
"""
constant = sp.symbols('c', real=True)
w_field_new = dg.Field(vector_field=sp.Matrix([0, constant]),
                       space=space)
print(w_field_new.is_killing()) # True

# 2.2 Find all functions so that the Lie derivative is equal to the metric
lie_derivative = w_field.lie_derivative()
print(lie_derivative) 
"""
For this to be equal to the metric g exactlly, we need to solve the two equations:
(1) : 4 * f'(x1) = 2
(2) : 2*x1**2 * h'(x2) + 2*x1*(fx1) = x1**2

From observation we that first we satisfy (1) by: f'(x1) = 1/2
Then we have f(x1) = 1/2 x1 + c1
Plugging this into (2) we get
2*x1**2 * h'(x2) + 2*x1*(1/2 x1) = x1**2 
2*x1**2 * h'(x2) + x1**2 + 2*x1*c1 = x1**2
Which we satisfy by c1=0 and h'(x2) = 0.
Then:
f(x1) = 1/2 x1
h(x2) = constant
We test this with the following
"""
f_test = sp.Rational(1,2) * x1
h_test = constant
W_field_lie = dg.Field(vector_field=sp.Matrix([f_test, h_test]),
                       space=space)
print(W_field_lie.lie_derivative())
print(W_field_lie.lie_derivative() == space.metric) # This returns true
"""
Thus we have
w_field = [1/2 x1, c] for any real number c. 
"""





