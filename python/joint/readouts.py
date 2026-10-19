"""Exact readouts in the notation of the current manuscript.

Paper correspondence (Section 6):
    normalization(n) = C_n                    [Lemma 6.9]
    B(n,k), Bprime(n,k) = B_{n,k}, B'_{n,k}    [Lemmas 6.9–6.11]
    forward_values(n) = (bold q, bold p3, bold p4)
    dual_rows(n) = (bold a, bold b, bold d)    [Lemma 6.12]
    Lam(n) = Lambda_n                        [Lemma 6.13]

All returned sums use integers or fractions.Fraction. Bprime is the
Gamma-interpolation derivative, including simple-zero tails; it is not
the derivative of a piecewise-zero definition. Notebook 2 supplies the
symbolic telescoping and boundary arguments for the all-index identities.
"""
from fractions import Fraction as Fr
from math import comb, factorial as fac, floor, ceil
from collections import Counter
import math
import mpmath as mp
import sympy

def vp(r, p):
    # p-adic valuation of a nonzero rational
    r = Fr(r)
    if r == 0:
        return math.inf
    a, b, v = abs(r.numerator), r.denominator, 0
    while a % p == 0:
        a //= p; v += 1
    while b % p == 0:
        b //= p; v -= 1
    return v

def legendre_bound(terms, p):
    # lower bound for v_p(sum of terms): the smallest valuation of a nonzero term
    return min(vp(t, p) for t in terms if t != 0)

_Hc = {}
def H(m, r=1):
    # H^{(r)}_m = sum_{i=1}^m i^{-r}, exact; H_0 = 0
    if (m, r) not in _Hc:
        _Hc[(m, r)] = sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))
    return _Hc[(m, r)]

def binom(a, b):
    return comb(a, b) if 0 <= b <= a else 0

def B(n, k):
    """Integer binomial summand B_{n,k}; ordinary support 0 <= k <= 14n."""
    return (binom(k + 12*n, 10*n) * binom(k + 13*n, 10*n) * binom(k + 14*n, 10*n)
            * binom(15*n, n + k) * binom(15*n, k))

def Bprime(n, k):
    """Exact derivative B'_{n,k} on -n <= k <= 15n, tails included."""
    if 0 <= k <= 14*n:
        br = (sum(H(k + (10 + u)*n) - H(k + u*n) for u in (2, 3, 4))
              + H(14*n - k) - H(n + k) + H(15*n - k) - H(k))
        return B(n, k) * br
    if 14*n < k <= 15*n:
        num = binom(k + 12*n, 10*n) * binom(k + 13*n, 10*n) * binom(k + 14*n, 10*n) * binom(15*n, k)
        return Fr((-1)**k * num, (n + k) * binom(n + k - 1, 15*n))
    if -n <= k < 0:
        i = -k
        num = binom(12*n - i, 10*n) * binom(13*n - i, 10*n) * binom(14*n - i, 10*n) * binom(15*n, n - i)
        return Fr((-1)**(i - 1) * num, (15*n + i) * binom(15*n + i - 1, i - 1))
    return Fr(0)

def forward_terms(n):
    # the terms of the three sums, as exact rationals
    q = [Fr(B(n, k)) for k in range(0, 14*n + 1)]
    p = {s: [Fr(B(n, k)) * H(4*n + k, s) for k in range(-n, 15*n + 1)]
            + [-Fr(1, s - 1) * Bprime(n, k) * H(4*n + k, s - 1) for k in range(-n, 15*n + 1)]
         for s in (3, 4)}
    return {"q": q, "p3": p[3], "p4": p[4]}

def Bstar(n, j):
    """Auxiliary triple-pole binomial product for the dual kernel."""
    return (binom(16*n + j, 14*n) * binom(17*n + j, 16*n) * binom(11*n, j)
            * binom(11*n, j - n) * binom(8*n + 1, j - 2*n))

def Nstar(n):
    """Auxiliary factorial normalizer for the dual leading Laurent term."""
    return Fr(fac(11*n)**2 * fac(8*n + 1), fac(14*n) * fac(16*n))

def Bstar_prime(n, j):
    """Gamma derivative of the dual binomial product, holding (-1)^j fixed."""
    if 2*n <= j <= 10*n + 1:
        br = (H(16*n + j) - H(2*n + j) + H(17*n + j) - H(n + j) + H(11*n - j) - H(j)
              + H(12*n - j) - H(j - n) + H(10*n + 1 - j) - H(j - 2*n))
        return Bstar(n, j) * br
    if n <= j <= 11*n:
        P = binom(16*n + j, 14*n) * binom(17*n + j, 16*n) * binom(11*n, j) * binom(11*n, j - n)
        if j < 2*n:
            tail = Fr(fac(8*n + 1) * fac(2*n - j - 1), fac(10*n + 1 - j))
        else:
            tail = Fr(fac(8*n + 1) * fac(j - 10*n - 2), fac(j - 2*n))
        return (-1)**(j + 1) * P * tail
    return Fr(0)

def laurent(n):
    # Laurent coefficients rho[j][i] = coefficient of eps^{-i} in R_n(-j+eps), at every pole s = -j
    num = [1 + 2*n + i for i in range(14*n)] + [1 + n + i for i in range(16*n)]
    den = Counter([-i for i in range(11*n + 1)] + [-n - i for i in range(11*n + 1)]
                  + [-2*n - i for i in range(8*n + 2)])
    out = {}
    for rho, m in den.items():
        j = -rho
        C, L1, L2 = Fr((-1)**n), Fr(0), Fr(0)
        for r in num:                     # factor (r - s) = (r + j)(1 - eps/(r + j))
            dd = Fr(r + j); C *= dd; L1 -= 1/dd; L2 -= 1/(2*dd*dd)
        for r2, e in den.items():         # factor (s - r2) = (-j - r2)(1 + eps/(-j - r2))
            if r2 == rho:
                continue
            dd = Fr(-j - r2); C /= dd**e; L1 -= e/dd; L2 += e/(2*dd*dd)
        g = [C, C*L1, C*(L1*L1/2 + L2)]
        out[j] = (m, {m: g[0], m - 1: g[1], m - 2: g[2]})
    return out

def Lam(n):
    """Paper Lambda_n; includes 1/6 and no residual endpoint Delta factor."""
    f = fac
    return Fr(f(8*n) * f(9*n)**2 * f(11*n)**3 * f(12*n)**3 * f(13*n)**2,
              6 * f(16*n)**2 * f(17*n)**2 * f(18*n)**2 * f(19*n))

def dual_terms(n):
    Ns = Nstar(n)
    a = [3 * (-1)**j * Fr(j)**k * Bstar(n, j) / Ns for k in range(3) for j in range(2*n, 10*n + 2)]
    b = ([-2 * (-1)**j * (k * Fr(j)**(k - 1) if k > 0 else 0) * Bstar(n, j) / Ns
          for k in range(3) for j in range(n, 11*n + 1)]
         + [-2 * (-1)**j * Fr(j)**k * Bstar_prime(n, j) / Ns for k in range(3) for j in range(n, 11*n + 1)])
    d = []
    for j, (m, rho) in laurent(n).items():
        if m == 3:   # sanity: the top coefficient at a triple pole is (-1)^j B*/N*
            assert rho[3] == (-1)**j * Bstar(n, j) / Ns
        for k in range(3):
            for i in range(1, m + 1):
                c = sum(comb(k, l) * Fr(j)**(k - l) * (-1)**l * rho.get(i + l, 0) for l in range(k + 1))
                d.append(-i * c * H(2*n + j, i + 1))
    return {"a": a, "b": b, "d": d}

def normalization(n):
    """Paper C_n; bold readouts equal C_n/12 times matrix readouts."""
    return binom(12*n,10*n)*binom(13*n,10*n)*binom(14*n,10*n)*binom(15*n,n)

def forward_values(n):
    """Return [bold q_n, bold p_{3,n}, bold p_{4,n}] as exact Fractions."""
    terms = forward_terms(n)
    return [sum(terms[k], Fr(0)) for k in ("q", "p3", "p4")]

def dual_rows(n):
    """Return [bold a_n, bold b_n, bold d_n], including the n=0 simple pole."""
    result = [[Fr(0) for k in range(3)] for r in range(3)]
    for j, (order, rho) in laurent(n).items():
        for k in range(3):
            c = {i: sum((comb(k,l)*Fr(j)**(k-l)*(-1)**l*rho.get(i+l,0)
                         for l in range(k+1)),Fr(0)) for i in range(1,4)}
            result[0][k] += 3*c[3]
            result[1][k] += 2*c[2]
            result[2][k] -= sum(i*c[i]*H(2*n+j,i+1) for i in range(1,4))
    return result

def lcm_to(m):
    """Paper d_m = lcm(1,...,m), with d_0=1."""
    return math.lcm(*range(1,m+1)) if m else 1
