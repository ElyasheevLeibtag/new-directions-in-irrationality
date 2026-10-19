"""Exact prescribed creative telescoping, using FLINT over Q[n][k].

No recurrence fitting, interpolation, Maple, or numerical arithmetic is used.
The notebook explains the equations and the one resonant coefficient.
"""
from __future__ import annotations

import ast
import json
import math
import time
from pathlib import Path

from flint import fmpq, fmpq_poly

N = fmpq_poly([0, 1])
ZERO = fmpq_poly([])
ONE = fmpq_poly([1])


class RationalPolynomial:
    """Small exact Q(n) class, also used by the restricted JSON parser."""
    def __init__(self, numerator=0, denominator=1):
        if isinstance(numerator, RationalPolynomial):
            self.p, self.q = numerator.p, numerator.q
            return
        self.p = numerator if isinstance(numerator, fmpq_poly) else fmpq_poly([numerator])
        self.q = denominator if isinstance(denominator, fmpq_poly) else fmpq_poly([denominator])
        if not self.q:
            raise ZeroDivisionError("zero rational-function denominator")
        g = self.p.gcd(self.q)
        self.p, self.q = self.p // g, self.q // g
        c = self.q.leading_coefficient()
        self.p, self.q = self.p / c, self.q / c

    def __add__(self, other):
        other = RationalPolynomial(other)
        if self.q == other.q:
            return RationalPolynomial(self.p + other.p, self.q)
        return RationalPolynomial(self.p * other.q + other.p * self.q, self.q * other.q)

    __radd__ = __add__

    def __neg__(self):
        return RationalPolynomial(-self.p, self.q)

    def __sub__(self, other):
        return self + (-RationalPolynomial(other))

    def __rsub__(self, other):
        return RationalPolynomial(other) - self

    def __mul__(self, other):
        other = RationalPolynomial(other)
        return RationalPolynomial(self.p * other.p, self.q * other.q)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = RationalPolynomial(other)
        return RationalPolynomial(self.p * other.q, self.q * other.p)

    def __rtruediv__(self, other):
        return RationalPolynomial(other) / self

    def __pow__(self, power):
        if power < 0:
            return RationalPolynomial(self.q ** (-power), self.p ** (-power))
        return RationalPolynomial(self.p ** power, self.q ** power)

    def shift(self, h):
        return RationalPolynomial(self.p(N + h), self.q(N + h))


def parse_rational(expression):
    """Accept integers, n, +, -, *, /, and integer powers; no eval()."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return RationalPolynomial(node.value)
        if isinstance(node, ast.Name) and node.id == "n":
            return RationalPolynomial(N)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = visit(node.operand)
            return value if isinstance(node.op, ast.UAdd) else -value
        if isinstance(node, ast.BinOp):
            left = visit(node.left)
            if isinstance(node.op, ast.Pow):
                exponent = node.right
                sign = 1
                if isinstance(exponent, ast.UnaryOp) and isinstance(exponent.op, ast.USub):
                    exponent, sign = exponent.operand, -1
                if not isinstance(exponent, ast.Constant) or type(exponent.value) is not int:
                    raise ValueError("only integer powers are allowed")
                return left ** (sign * exponent.value)
            right = visit(node.right)
            if isinstance(node.op, ast.Add): return left + right
            if isinstance(node.op, ast.Sub): return left - right
            if isinstance(node.op, ast.Mult): return left * right
            if isinstance(node.op, ast.Div): return left / right
        raise ValueError("unsupported expression in matrix JSON")
    return visit(ast.parse(str(expression), mode="eval").body)


def polynomial_lcm(a, b):
    return (a // a.gcd(b)) * b


def primitive_coefficients(values):
    """Return primitive integer polynomials in the supplied shift order."""
    denominator = ONE
    for value in values:
        denominator = polynomial_lcm(denominator, value.q)
    coefficients = [value.p * (denominator // value.q) for value in values]
    common = coefficients[0]
    for value in coefficients[1:]:
        common = common.gcd(value)
    coefficients = [value // common for value in coefficients]
    integer_denominator = math.lcm(*(int(c.denom()) for c in coefficients))
    coefficients = [value * integer_denominator for value in coefficients]
    content = math.gcd(*(int(c) for value in coefficients for c in value.coeffs()))
    coefficients = [value / content for value in coefficients]
    if coefficients[0].leading_coefficient() < 0:
        coefficients = [-value for value in coefficients]
    return coefficients


def factorial_step(a):
    result = ONE
    for j in range(1, a + 1):
        result *= a * N + j
    return result


def prescribed_recurrence(data_directory, kind):
    """Descending shifts: B[0]s(n+3)+...+B[3]s(n)=0."""
    data_directory = Path(data_directory)
    file = "matrices_for_zeta3_4.json" if kind == "forward" else "dual_companion.json"
    C = json.loads((data_directory / file).read_text())["C"]
    g = [parse_rational(C[j][2]) for j in range(3)]
    if kind == "forward":
        rho = RationalPolynomial(
            factorial_step(12) * factorial_step(13) * factorial_step(15),
            factorial_step(10) ** 3 * factorial_step(2) * factorial_step(3)
            * factorial_step(4) * factorial_step(1),
        )
        for j in range(3):
            for h in range(j, 3):
                g[j] *= rho.shift(h)
    elif kind != "dual":
        raise ValueError("kind must be forward or dual")
    coefficients = primitive_coefficients([RationalPolynomial(1), -g[2], -g[1], -g[0]])
    for j in range(3):
        assert not (RationalPolynomial(coefficients[3 - j], coefficients[0]) + g[j]).p
    return coefficients


def trim(p):
    while p and not p[-1]:
        p.pop()
    return p


def add(p, q):
    result = [ZERO] * max(len(p), len(q))
    for j, value in enumerate(p): result[j] += value
    for j, value in enumerate(q): result[j] += value
    return trim(result)


def scale(p, c):
    return trim([value * c for value in p])


def multiply(p, q):
    if not p or not q: return []
    result = [ZERO] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i + j] += a * b
    return trim(result)


def linear_product(factors):
    """Factors are (coefficient of n, coefficient of k, constant)."""
    result = [ONE]
    for a, b, c in factors:
        result = multiply(result, [a * N + c, fmpq_poly([b])])
    return result


def shifted_k(p):
    result = [ZERO] * len(p)
    for i, coefficient in enumerate(p):
        for j in range(i + 1):
            result[j] += coefficient * math.comb(i, j)
    return trim(result)


def apply_operator(A, B, P):
    return add(multiply(A, shifted_k(P)), scale(multiply(B, P), -1))


def solve_prescribed(A, B, F, top_degree, progress=None):
    """Solve A(k)P(k+1)-B(k)P(k)=F by descending exact coefficients.

    The diagonal pivot for k^i is i-top_degree. All pivots but the
    leading one are nonzero integers. The remaining low equations fix
    that one coefficient by exact polynomial division in Q[n].
    """
    d = top_degree
    assert len(A) == len(B) == 6 and A[5] == B[5] == ONE
    assert A[4] - B[4] == fmpq_poly([-d])
    assert len(F) <= d + 4
    residual = F[:] + [ZERO] * (d + 5 - len(F))
    basis = [math.comb(d, j) * ONE for j in range(d + 1)]
    homogeneous_residual = add(multiply(A, basis), [ZERO] * d + scale(B, -1))
    homogeneous_residual += [ZERO] * (d + 5 - len(homogeneous_residual))
    particular = [ZERO] * (d + 1)
    homogeneous = [ZERO] * d + [ONE]
    for i in range(d - 1, -1, -1):
        pivot = i - d
        part = residual[i + 4] / pivot
        hom = -homogeneous_residual[i + 4] / pivot
        particular[i], homogeneous[i] = part, hom
        for j in range(i + 5):
            entry = ZERO
            for a in range(6):
                if 0 <= j - a <= i:
                    entry += A[a] * math.comb(i, j - a)
            if 0 <= j - i < 6:
                entry -= B[j - i]
            residual[j] -= part * entry
            homogeneous_residual[j] += hom * entry
        assert not residual[i + 4] and not homogeneous_residual[i + 4]
        if progress and i % 20 == 0:
            progress(f"Solved coefficients through k^{i}")
    index = next((j for j in range(4) if homogeneous_residual[j]), None)
    if index is None:
        raise ArithmeticError("remaining equation does not determine the leading coefficient")
    top, remainder = divmod(residual[index], homogeneous_residual[index])
    if remainder:
        raise ArithmeticError("certificate is not polynomial in n")
    for j in range(4):
        assert residual[j] == top * homogeneous_residual[j]
    P = [particular[j] + top * homogeneous[j] for j in range(d + 1)]
    # Independent full substitution: do not trust the triangular solver.
    assert apply_operator(A, B, P) == trim(F[:]), "exact certificate residual is nonzero"
    return trim(P)


def forward_problem(B):
    A = linear_product([(12, 1, 1), (13, 1, 1), (14, 1, 1), (14, -1, 42), (15, -1, 45)])
    Bminus = linear_product([(0, 1, 0), (1, 1, 3), (2, 1, 6), (3, 1, 9), (4, 1, 12)])
    D = linear_product([(a, 1, r) for a in range(1, 5) for r in range(1, 3 * a + 1)]
                       + [(14, -1, r) for r in range(1, 43)]
                       + [(15, -1, r) for r in range(1, 46)])
    K = ONE
    for r in range(1, 31): K *= (10 * N + r) ** 3
    F = []
    for j in range(4):
        normalization = ONE
        for r in range(1, 15 * j + 1): normalization *= (15 * N + r) ** 2
        for r in range(10 * j + 1, 31): normalization *= (10 * N + r) ** 3
        factors = ([(a, 1, r) for a in range(12, 15) for r in range(1, a * j + 1)]
                   + [(a, 1, r) for a in range(1, 5) for r in range(a * j + 1, 3 * a + 1)]
                   + [(14, -1, r) for r in range(14 * j + 1, 43)]
                   + [(15, -1, r) for r in range(15 * j + 1, 46)])
        F = add(F, scale(linear_product(factors), B[3 - j] * normalization))
    return A, Bminus, F, D, K


def dual_problem(h):
    D = linear_product([(a, 1, r) for a in range(1, 3) for r in range(1, h * a + 1)]
                       + [(b, -1, r) for b in range(11, 13) for r in range(1, h * b + 1)]
                       + [(10, -1, 1 + r) for r in range(1, 10 * h + 1)])
    A = scale(linear_product([(16, 1, 1), (17, 1, 1), (11, -1, 11 * h),
                             (12, -1, 12 * h), (10, -1, 1 + 10 * h)]), -1)
    Bminus = linear_product([(0, 1, 0), (-1, 1, 0), (-2, 1, 0), (1, 1, h), (2, 1, 2 * h)])
    return A, Bminus, D


def dual_shift(j, h):
    factors = ([(a, 1, r) for a in range(16, 18) for r in range(1, a * j + 1)]
               + [(a, 1, r) for a in range(1, 3) for r in range(a * j + 1, h * a + 1)]
               + [(b, -1, r) for b in range(11, 13) for r in range(b * j + 1, h * b + 1)]
               + [(10, -1, 1 + r) for r in range(10 * j + 1, 10 * h + 1)]
               + [(-a, 1, -r) for a in range(1, 3) for r in range(a * j)])
    return linear_product(factors)


def degree_pair(P):
    return max(c.degree() for c in P), len(P) - 1


def no_integer_roots_at_or_above(poly, start=0):
    """Prove nonvanishing at all integers >=start by a shifted sign test.

    The translated constant term must be nonzero; otherwise a root at
    start+1 would pass the nonnegative-coefficient test.
    """
    if not poly: return False
    if poly.leading_coefficient() < 0: poly = -poly
    return (poly(start) != 0 and poly(start + 1) != 0
            and all(c >= 0 for c in poly(N + start + 1).coeffs()))


def evaluate_k(P, k):
    result = ZERO
    for coefficient in reversed(P):
        result = result * k + coefficient
    return result


def generate_all(data_directory="data", progress=print):
    """Generate five certificates; every return has passed full substitution."""
    result = {}
    started = time.monotonic()
    B = prescribed_recurrence(data_directory, "forward")
    if progress: progress(f"Forward recurrence degrees: {[p.degree() for p in B]}")
    A, Bminus, F, D, K = forward_problem(B)
    P = solve_prescribed(A, Bminus, F, 114, progress)
    assert no_integer_roots_at_or_above(B[0])
    assert no_integer_roots_at_or_above(evaluate_k(D, 1 - 12 * N), 2)
    assert 5 + len(P) - 1 - (len(D) - 1) - 2 <= 0
    result["forward"] = dict(B=B, P=P, degrees=degree_pair(P))
    if progress: progress(f"Forward certificate verified; degrees {degree_pair(P)}")
    DB = prescribed_recurrence(data_directory, "dual")
    A, Bminus, D = dual_problem(3)
    F = []
    for j in range(4): F = add(F, scale(dual_shift(j, 3), DB[3 - j]))
    P = solve_prescribed(A, Bminus, F, 107, progress)
    assert no_integer_roots_at_or_above(DB[0])
    assert 5 + len(P) - 1 - (len(D) - 1) - 4 <= 0
    result["dual_scalar"] = dict(B=DB, P=P, degrees=degree_pair(P))
    if progress: progress(f"Dual scalar certificate verified; degrees {degree_pair(P)}")
    Qraw = json.loads((Path(data_directory) / "dual_matrices_zeta3_4.json").read_text())["G"]
    Q = [[parse_rational(Qraw[i][j]) for j in range(3)] for i in range(3)]
    A, Bminus, D = dual_problem(1)
    shifted = dual_shift(1, 1)
    for j in range(3):
        denominator = ONE
        for i in range(3): denominator = polynomial_lcm(denominator, Q[i][j].q)
        F = [ZERO] * j + scale(shifted, denominator)
        for i in range(3):
            coefficient = Q[i][j].p * (denominator // Q[i][j].q)
            F = add(F, [ZERO] * i + scale(D, -coefficient))
        P = solve_prescribed(A, Bminus, F, 35, progress)
        assert no_integer_roots_at_or_above(denominator)
        assert no_integer_roots_at_or_above(evaluate_k(D, 1 - 16 * N), 1)
        assert 5 + len(P) - 1 - (len(D) - 1) - 4 <= 0
        result[f"dual_row_{j}"] = dict(P=P, denominator=denominator, degrees=degree_pair(P))
        if progress: progress(f"Dual row {j} certificate verified; degrees {degree_pair(P)}")
    result["elapsed_seconds"] = round(time.monotonic() - started, 2)
    return result


if __name__ == "__main__":
    output = generate_all(Path(__file__).parent / "data")
    print({name: row["degrees"] for name, row in output.items() if isinstance(row, dict)})
    print("All five exact substitutions passed in", output["elapsed_seconds"], "seconds")
