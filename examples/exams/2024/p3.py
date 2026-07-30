import diffgeom as dg
import sympy as sp

# Exam 2024 Problem 1

def taylor(f:sp.Matrix, var:sp.symbols, n):
    func = 0
    for term in range(n):
        
        func += (f.diff(var))
