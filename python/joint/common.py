"""Small exact-arithmetic helpers shared by the proof notebooks."""
from pathlib import Path
from functools import lru_cache
import ast
import json
import sys
import sympy as sp
from ramanujantools import Matrix
from ramanujantools.flint_core import SymbolicMatrix, flint_ctx

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parent
n = sp.Symbol("n")
CTX = flint_ctx({n}, fmpz=True)

def expression(text):
    """Read rational functions in n; reject calls and other Python syntax."""
    def visit(t):
        if isinstance(t, ast.Constant) and type(t.value) is int:
            return sp.Integer(t.value)
        if isinstance(t, ast.Name) and t.id == "n":
            return n
        if isinstance(t, ast.UnaryOp) and isinstance(t.op, (ast.USub, ast.UAdd)):
            v = visit(t.operand)
            return -v if isinstance(t.op, ast.USub) else v
        if isinstance(t, ast.BinOp):
            a, b = visit(t.left), visit(t.right)
            if isinstance(t.op, ast.Add): return a + b
            if isinstance(t.op, ast.Sub): return a - b
            if isinstance(t.op, ast.Mult): return a * b
            if isinstance(t.op, ast.Div): return a / b
            if isinstance(t.op, ast.Pow) and b.is_Integer: return a ** b
        raise ValueError("Expected a rational expression in n")
    return visit(ast.parse(str(text), mode="eval").body)

@lru_cache(None)
def data(filename):
    return json.loads((ROOT / "data" / filename).read_text())

def matrix(filename, key):
    return Matrix([[expression(v) for v in row] for row in data(filename)[key]])

def sm(A):
    return SymbolicMatrix.from_sympy(Matrix(A), CTX)

def equal(A, B):
    return sm(A) == sm(B)

def scalar_equal(a, b):
    return equal(Matrix([[a]]), Matrix([[b]]))

def determinant(A):
    """Exact small determinant with FLINT rational-function arithmetic."""
    from itertools import permutations
    S = sm(A)
    d = A.rows
    total = sm(Matrix([[0]]))[0, 0]
    for p in permutations(range(d)):
        term = sm(Matrix([[1]]))[0, 0]
        for i, j in enumerate(p): term *= S[i, j]
        inversions = sum(p[i] > p[j] for i in range(d) for j in range(i+1, d))
        total += (-1 if inversions % 2 else 1) * term
    return SymbolicMatrix(1, 1, [total], CTX).factor()[0, 0]

def regular_nonnegative(expr):
    """Prove p(k)!=0 for every integer k>=0 by a shifted sign test.

    Check k=0 separately. The shifted polynomial has coefficients of one
    sign and a nonzero constant term, so it never vanishes for real n>=0.
    The constant-term check is essential: p(n)=n-1 must be rejected.
    """
    P = sp.Poly(expr, n)
    c = sp.Poly(P.as_expr().subs(n, n+1), n).all_coeffs()
    return P.eval(0) != 0 and c[-1] != 0 and (all(x >= 0 for x in c) or all(x <= 0 for x in c))

def step(a):
    return sp.prod(a*n+j for j in range(1, a+1))

def factorial_ratio_step(exponents):
    return sp.prod(step(a)**e for a,e in exponents.items())
