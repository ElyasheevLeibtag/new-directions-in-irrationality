"""Readable exact helpers for the paper's ζ(2) trajectory.

The public interface follows the PAPER upper-parameter order
(x_0,x_1,x_2,x_3) = (0,1,1,1) at n=0.  No Lean is imported or generated.
The independent archival certificates retain their original parameter order.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, factorial, isqrt, lcm, prod
from pathlib import Path
import json
import subprocess
import sys
import time

from flint import fmpq, fmpq_mat, fmpq_poly

START = (0, 1, 1, 1, 1, 1, 1)
DIRECTION = (-11, 10, 9, 11, 2, 1, 3)
IDENTITY = fmpq_mat([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
READOUT_ROWS = fmpq_mat([[0, 2, -1], [2, 0, 0]])  # p row, then q row
POLY_S = fmpq_poly([0, 1])
SHELLS = tuple((Fraction(a, b), Fraction(c, d)) for a, b, c, d in (
    (13, 2, 8, 1), (13, 3, 9, 2), (13, 4, 10, 3),
    (13, 5, 8, 3), (13, 6, 9, 4), (13, 12, 9, 8)))


def endpoint(n):
    """s_0+n t_0, in the paper's parameter order."""
    return tuple(s + n*t for s, t in zip(START, DIRECTION))


def theta(vertex):
    """Matrix of θ on row states (F,θF,θ²F), at z=1."""
    upper = vertex[:4]
    lower_minus_one = [y-1 for y in vertex[4:]]
    elementary = lambda values, k: sum(prod(t) for t in combinations(values, k))
    coefficients = (
        -prod(upper),
        prod(lower_minus_one) - elementary(upper, 3),
        elementary(lower_minus_one, 2) - elementary(upper, 2),
        sum(lower_minus_one) - sum(upper),
    )
    q0, q1, q2, q3 = coefficients
    if not q3:
        raise ZeroDivisionError(f"Singular differential-equation vertex: {vertex}")
    return fmpq_mat([[0, 0, fmpq(-q0, q3)],
                     [1, 0, fmpq(-q1, q3)],
                     [0, 1, fmpq(-q2, q3)]])


def unit_move(vertex, axis, sign):
    """Return the endpoint and exact contiguity matrix of one lattice move.

    Upper parameters can increase or decrease. Only lower increases are used
    by these routes; reject other lower moves instead of silently misusing a
    formula. Matrices multiply on the right of row states.
    """
    if sign not in (-1, 1) or (axis >= 4 and sign != 1):
        raise ValueError("This helper supports upper ±1 and lower +1 moves.")
    target = list(vertex)
    target[axis] += sign
    target = tuple(target)
    if axis < 4 and sign == 1:
        matrix = IDENTITY + theta(vertex)*fmpq(1, vertex[axis])
    elif axis < 4:
        matrix = (IDENTITY + theta(target)*fmpq(1, target[axis])).inv()
    else:
        matrix = (IDENTITY + theta(target)*fmpq(1, vertex[axis])).inv()
    return target, matrix


def four_block_walk(n, plateau=3, backend="native"):
    """Compute Φ_M(n) by the A,B,C,D route from Appendix app:z2-blocks.

    plateau=3 is the main-text path; plateau=4 is an independent path check.
    backend='ramanujantools' constructs every unit move with pFq(4,3,1).
    All entries and comparisons are rational and exact.
    """
    if n < 0 or plateau not in (3, 4):
        raise ValueError("Use n≥0 and plateau in {3,4}.")
    if backend == "ramanujantools":
        import sympy as sp
        from ramanujantools import Position
        from ramanujantools.cmf.pfq import pFq
        axes = sp.symbols("x0 x1 x2 x3 y0 y1 y2")
        cmf = pFq(4, 3, 1)
    elif backend != "native":
        raise ValueError("backend must be 'native' or 'ramanujantools'.")
    blocks = (
        (plateau*n, (0, 1, 1, 1, 0, 0, 0)),
        (n, (0, 0, 0, 0, 2, 1, 3)),
        (n, (0, 10-plateau, 9-plateau, 11-plateau, 0, 0, 0)),
        (11*n, (-1, 0, 0, 0, 0, 0, 0)),
    )
    vertex, walk = START, IDENTITY
    for count, direction in blocks:
        for _ in range(count):
            for axis, amount in enumerate(direction):
                sign = 1 if amount >= 0 else -1
                for _ in range(abs(amount)):
                    if backend == "native":
                        vertex, matrix = unit_move(vertex, axis, sign)
                    else:
                        substitutions = Position(dict(zip(axes, vertex)))
                        matrix = cmf.M(axes[axis], sign > 0).subs(substitutions)
                        matrix = fmpq_mat([[fmpq(str(matrix[i, j])) for j in range(3)]
                                           for i in range(3)])
                        changed = list(vertex)
                        changed[axis] += sign
                        vertex = tuple(changed)
                    walk = walk*matrix
    assert vertex == endpoint(n)
    return walk


def matrix_readout(n):
    """Return (p_n,q_n) with ρ_p=(0,2,-1), ρ_q=(2,0,0)."""
    readout = READOUT_ROWS*four_block_walk(n)
    return readout[0, 0], readout[1, 0]


def normalization_C(n):
    """C_n=∏_{u=9,10,11} binomial(un,8n)."""
    return prod(comb(u*n, 8*n) for u in (9, 10, 11))


def binomial_polynomial(n):
    """The polynomial F(s)=∏ binomial(s+un,8n), not the function F_n."""
    result = fmpq_poly([1])
    for u in (9, 10, 11):
        for k in range(8*n):
            result *= POLY_S + u*n-k
        result /= factorial(8*n)
    return result


def divided_sum(n, polynomial=None):
    """S_0(s); exact division fills the removable singularity at each j."""
    polynomial = binomial_polynomial(n) if polynomial is None else polynomial
    result = fmpq_poly([])
    for j in range(11*n+1):
        quotient, remainder = divmod(polynomial-polynomial(j), POLY_S-j)
        assert not remainder
        result += (-1)**j*comb(11*n, j)*quotient
    return result


def J_eta(polynomial, eta):
    """J_η(P)=Σ_h (-1)^h Δ^hP(-η)/(h+1), a finite Newton-basis sum."""
    if not polynomial:
        return fmpq(0)
    differences = [polynomial(k-eta) for k in range(polynomial.degree()+1)]
    total = fmpq(0)
    for h in range(len(differences)):
        total += fmpq((-1)**h, h+1)*differences[0]
        differences = [b-a for a, b in zip(differences, differences[1:])]
    return total


def closed_form_readout(n, eta=None):
    """Return (C_n p_n,C_n q_n), as printed in Appendix app:z2-blocks.

    n=0 uses η=0. For n≥1, any integer 2n≤η≤9n is permitted.
    """
    eta = 2*n if eta is None else eta
    if not 2*n <= eta <= 9*n:
        raise ValueError("η must be an integer between 2n and 9n.")
    polynomial = binomial_polynomial(n)
    coefficients = [(-1)**j*comb(11*n, j)*polynomial(j) for j in range(11*n+1)]
    harmonics = [fmpq(0)]
    for k in range(1, 11*n+eta+1):
        harmonics.append(harmonics[-1] + fmpq(1, k*k))
    cp = 2*(sum((c*harmonics[j+eta] for j, c in enumerate(coefficients)), fmpq(0))
            + J_eta(divided_sum(n, polynomial), eta))
    cq = 2*sum(coefficients, fmpq(0))
    return cp, cq


def lcm_to(m):
    """d_m=lcm(1,...,m), with d_0=1."""
    return lcm(*range(1, m+1)) if m else 1


def primes_to(m):
    """Simple sieve, adequate for the explicit finite checks."""
    sieve = [True]*(m+1)
    for p in range(2, isqrt(m)+1):
        if sieve[p]:
            for composite in range(p*p, m+1, p):
                sieve[composite] = False
    return [p for p in range(2, m+1) if sieve[p]]


def common_divisors(n):
    """Return (D_n^car,D_n^sh) with the paper's strict shell endpoints."""
    sys.path.insert(0, str(Path(__file__).parent / "certificates"))
    from rate_certificate import phi
    carry = shell = 1
    for p in primes_to(13*n):
        if p*p <= 50*n:
            continue
        carry *= p**phi(Fraction(n, p))
        if any(left*n < p < right*n for left, right in SHELLS):
            shell *= p
    return carry, shell


def run_certificate(relative_script):
    """Replay a supplied Python certificate in its own process; no cached pass.

    JSON certificates and logs are regenerated. The archived scripts write
    their JSON next to their source. The caller receives the parsed new JSON.
    """
    directory = Path(__file__).parent / "certificates"
    script = directory / relative_script
    output_names = {
        "analytic/certify_analytic_rates.py": "analytic/analytic_certificates.json",
        "check_finite_formula.py": "finite_formula_checks.json",
        "check_main_shell.py": "main_shell_checks.json",
        "check_fifth_shell.py": "fifth_shell_checks.json",
        "check_shell_independent.py": "independent_checks.json",
    }
    output = directory / output_names.get(relative_script, str(Path(relative_script).with_suffix(".json")))
    started = time.perf_counter()
    completed = subprocess.run([sys.executable, str(script)], cwd=directory,
                               text=True, capture_output=True)
    logs = Path(__file__).parent / "generated"
    logs.mkdir(exist_ok=True)
    (logs / (script.stem + ".log")).write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"{script.name} failed:\n{completed.stdout}\n{completed.stderr}")
    result = json.loads(output.read_text())
    return result, round(time.perf_counter()-started, 3)
