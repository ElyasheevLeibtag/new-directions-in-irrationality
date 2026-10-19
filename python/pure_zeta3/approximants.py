"""Exact partial fractions and approximants of Appendix app:zeta3-measure.

Q_n = sum_j a_j and P_{3,n} = sum_j(a_j H_{j+2n}^{(3)}
+ b_j H_{j+2n}^{(2)}/2). These capital symbols are the native Meijer-G
normalization, not the section's forward q_n and p_{3,n}.

The formulas below preserve the supplied partial_fractions.py arithmetic.
python-flint fmpq stores exact rational numbers, never binary floats.
"""
from flint import fmpq
import math

def partial_fractions_direct(n):
    """Compute a_j and b_j directly from kernel products; useful as an independent finite check."""
    numerator = [k for u in (0, 1, 2) for k in range(u * n + 1, (10 + u) * n + 1)]
    den = [*range(n, 17 * n + 1), *range(2 * n, 16 * n + 1)]
    result = {}
    for j in range(n, 17 * n + 1):
        multiplicity = den.count(j)
        c = fmpq(1)
        for k in numerator:
            c *= -j - k
        for h in den:
            if h != j:
                c /= h - j
        if multiplicity == 1:
            result[j] = (fmpq(0), c)
        else:
            L = sum((fmpq(1, -j - k) for k in numerator), fmpq(0)) - sum((fmpq(1, h - j) for h in den if h != j), fmpq(0))
            result[j] = (c, c * L)
    assert sum((b for (a, b) in result.values()), fmpq(0)) == 0
    return result

def partial_fractions(n):
    """Compute a_j and b_j from exact factorial and harmonic formulas, including both simple-pole tails."""
    if n == 0:
        return {0: (fmpq(1), fmpq(0))}
    f = [1]
    H = [fmpq(0)]
    for k in range(1, 29 * n + 1):
        f.append(f[-1] * k)
        H.append(H[-1] + fmpq(1, k))
    cs = {}
    for j in range(n, 17 * n + 1):
        numer = f[j + 10 * n] * f[j + 11 * n] * f[j + 12 * n]
        denom = f[j] * f[j + n] * f[j + 2 * n] * f[j - n] * f[17 * n - j]
        if j < 2 * n:
            a = fmpq(0)
            b = (-1) ** (j - n) * fmpq(numer * f[2 * n - j - 1], denom * f[16 * n - j])
        elif j > 16 * n:
            a = fmpq(0)
            b = (-1) ** (j - n + 1) * fmpq(numer * f[j - 16 * n - 1], denom * f[j - 2 * n])
        else:
            a = (-1) ** n * fmpq(numer, denom * f[j - 2 * n] * f[16 * n - j])
            L = sum((H[j + u * n] - H[j + (10 + u) * n] for u in (0, 1, 2)), fmpq(0)) + H[j - n] - H[17 * n - j] + H[j - 2 * n] - H[16 * n - j]
            b = a * L
        cs[j] = (a, b)
    assert sum((b for (a, b) in cs.values()), fmpq(0)) == 0
    return cs

def kernel_readout(n, eta=None):
    """Return (Q_n,P_{3,n},P_{4,n}) exactly, with the paper cutoff eta=2n by default."""
    if eta is None:
        eta = 2 * n
    cs = partial_fractions(n)
    H = {r: [fmpq(0)] for r in (2, 3, 4)}
    for k in range(1, 17 * n + eta + 1):
        for r in H:
            H[r].append(H[r][-1] + fmpq(1, k ** r))
    Q = sum((a for (a, b) in cs.values()), fmpq(0))
    P3 = sum((a * H[3][j + eta] + b * H[2][j + eta] / 2 for (j, (a, b)) in cs.items()), fmpq(0))
    P4 = sum((a * H[4][j + eta] + b * H[3][j + eta] / 3 for (j, (a, b)) in cs.items()), fmpq(0))
    return (Q, P3, P4)

def finite_integer_check(n, rows):
    """Finite regression of the uniform divisor, never its proof.

    The prime conditions are exactly those in Appendix B.4:
    p>29, p^2>30n, p does not divide n. Values on Farey walls are
    impossible for an eligible prime. For x>=1 use the periodic tail.
    """
    from fractions import Fraction as F
    from sympy import primerange
    from bisect import bisect_right
    if n < 1:
        raise ValueError('The divisor test uses n >= 1.')
    ends = [F(row['hi']) for row in rows]
    (divisor, first_count, second_count) = (1, 0, 0)
    for p in primerange(max(29, math.isqrt(30 * n)) + 1, 19 * n + 1):
        p = int(p)
        if n % p == 0:
            continue
        x = F(n, p)
        fractional_x = x - int(x)
        row = rows[bisect_right(ends, fractional_x)]
        assert F(row['lo']) < fractional_x < F(row['hi'])
        baseline = row['g3'] if x < 1 else row['tail']
        first = int(row['valid'])
        second = int(row['second']['valid'])
        first_count += first
        second_count += second
        divisor *= p ** (baseline + first + second)
    (Q, P3, _) = kernel_readout(n)
    Cprime = fmpq(math.factorial(15 * n) ** 2, math.factorial(10 * n) ** 3)
    d19 = math.lcm(*range(1, 19 * n + 1))
    normalization = 6 * Cprime * d19 ** 3
    d16 = math.lcm(*range(1, 16 * n + 1))
    for (a, b) in partial_fractions(n).values():
        assert (Cprime * a).denominator == 1
        assert (Cprime * d16 * b).denominator == 1
    (q_integer, p_integer) = (normalization * Q / divisor, normalization * P3 / divisor)
    assert q_integer.denominator == p_integer.denominator == 1
    assert int(q_integer.numerator) % int((P3 / Q).denominator) == 0
    return dict(
        n=n,
        first_order_primes=first_count,
        second_order_primes=second_count,
        integer_denominator_digits=len(str(abs(q_integer.numerator))),
        reduced_denominator_digits=len(str((P3 / Q).denominator)),
    )
