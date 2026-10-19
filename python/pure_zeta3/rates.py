"""Certified rational rate calculations for Appendix B.4 (app:zeta3-measure).

The matrices, characteristic polynomials, signs, logarithm enclosures, and
profile integration are verified here. The applicability of these constants
to the selected readout uses the paper's separate asymptotic arguments.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as sp

START = (1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
DIRECTION = (16, 17, -11, -12, -10, 2, 1, 0, -1, -2)


def native_trajectory():
    """The opposite-direction Meijer G trajectory in Appendix B.4."""
    from ramanujantools import Position
    from ramanujantools.cmf import MeijerG
    axes = sp.symbols('a0:5') + sp.symbols('b0:5')
    return MeijerG(5, 2, 5, 5, 1).trajectory_matrix(
        Position(dict(zip(axes, DIRECTION))),
        Position(dict(zip(axes, START))))


def matrix_at_infinity(matrix):
    """Limit of S_n^(-1) M(n) S_(n+1), S_n=diag(1,n^-1,n^-2,n^-3)."""
    n = sp.Symbol('n')
    limit = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            numerator, denominator = sp.cancel(matrix[i, j]).as_numer_denom()
            if numerator == 0:
                continue
            a, b = sp.Poly(numerator, n), sp.Poly(denominator, n)
            degree = a.degree() - b.degree()
            assert degree <= j - i
            if degree == j - i:
                limit[i, j] = a.LC() / b.LC()
    reference = json.loads((Path(__file__).parent/'data/spectral_reference.json').read_text())
    assert limit == sp.Matrix(reference['limiting_matrix'])
    assert limit.charpoly().all_coeffs() == list(map(sp.Rational, reference['characteristic_coefficients']))
    return limit


def _log_series(q, terms=40):
    """For 1 <= q <= 2, log(q)=2 atanh((q-1)/(q+1)); exact tail bound."""
    z = (q - 1) / (q + 1)
    lower = sum((2*z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    tail = 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lower, lower+tail


def log_enclosure(q):
    q, power = F(q), 0
    while q >= 2:
        q /= 2
        power += 1
    while q < 1:
        q *= 2
        power -= 1
    lo, hi = _log_series(q)
    lo2, hi2 = _log_series(F(2))
    return (lo+power*lo2, hi+power*hi2) if power >= 0 else (lo+power*hi2, hi+power*lo2)


def spectral_enclosures(limit=None):
    """Bracket the positive root; infer complex-pair modulus from Vieta."""
    reference = json.loads((Path(__file__).parent/'data/spectral_reference.json').read_text())
    coefficients = list(map(F, reference['characteristic_coefficients']))
    if limit is not None:
        assert list(map(F, map(str, limit.charpoly().all_coeffs()))) == coefficients
    # The quartic is (X+1)(X^3+aX^2+bX+c).
    a = coefficients[1]-1
    b = coefficients[2]-a
    c = coefficients[3]-b
    assert c == coefficients[4]
    discriminant = a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    assert discriminant < 0
    lower = F(138836415024165002899112)
    upper = lower+1
    evaluate = lambda x: ((x+a)*x+b)*x+c
    assert evaluate(lower) < 0 < evaluate(upper)
    assert c < 0 and -c/lower < 1
    log_lambda = log_enclosure(lower)[0], log_enclosure(upper)[1]
    log_rho = log_enclosure(-c/upper)[0]/2, log_enclosure(-c/lower)[1]/2
    g = log_lambda[0]-log_rho[1], log_lambda[1]-log_rho[0]
    cprime = tuple(30*v for v in log_enclosure(F(3, 2)))
    scaled_height = tuple(log_lambda[i]+cprime[i] for i in (0, 1))
    scaled_error = tuple(log_rho[i]+cprime[i] for i in (0, 1))
    assert g[0] > F('112.807539247454')
    assert scaled_height[1] < F('65.451536565866')
    assert scaled_error[1] < F('-47.356002681589')
    return dict(raw_Q_growth=log_lambda, raw_linear_error=log_rho,
                Cprime_growth=cprime, g=g, scaled_Q_growth=scaled_height,
                scaled_linear_error=scaled_error)


def density_enclosures(rows, periods=1000):
    """Sum the exact floor profile with certified infinite periodic tails.

    Below x=1 use g3; above x=1 use tail. First-order and second-order
    indicators are periodic, including the first period. Fixed-point sums
    round each term downward at 10^-18; tail logs round at 10^-24.
    """
    from chambers import merge_runs
    scale, tail_scale = 10**18, 10**24
    def floor_scaled(x, denominator):
        return F(x.numerator*denominator//x.denominator, denominator)
    def small_log(q):
        z = (q-1)/(q+1)
        lo = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(4)), F(0))
        err = 2*z**9/(9*(1-z*z))
        return floor_scaled(lo, tail_scale), -floor_scaled(-(lo+err), tail_scale)
    def tail(a, b):
        lo = small_log((periods+b)/(periods+a))[0]
        hi = small_log((periods-F(1, 2)+b)/(periods-F(1, 2)+a))[1]
        lo += (1/(periods+a)-1/(periods+b))/2
        return floor_scaled(lo, tail_scale), hi
    totals = {key: [F(0), F(0)] for key in ('baseline', 'first_order', 'second_order')}
    for row in rows:
        a, b = F(row['lo']), F(row['hi'])
        if row['g3']:
            value = row['g3']*(1/a-1/b)
            lo = floor_scaled(value, scale)
            totals['baseline'][0] += lo
            totals['baseline'][1] += lo+F(1, scale)
        if row['tail'] or row['valid']:
            tail_lo, tail_hi = tail(a, b)
        if row['tail']:
            lo = sum((floor_scaled(row['tail']*(1/(t+a)-1/(t+b)), scale)
                      for t in range(1, periods)), F(0))
            totals['baseline'][0] += lo+row['tail']*tail_lo
            totals['baseline'][1] += lo+F(periods-1, scale)+row['tail']*tail_hi
        if row['valid']:
            lo = sum((floor_scaled(1/(t+a)-1/(t+b), scale)
                      for t in range(periods)), F(0))
            totals['first_order'][0] += lo+tail_lo
            totals['first_order'][1] += lo+F(periods, scale)+tail_hi
    for a, b in merge_runs(rows, second=True):
        lo = sum((floor_scaled(1/(t+a)-1/(t+b), scale)
                  for t in range(periods)), F(0))
        tail_lo, tail_hi = tail(a, b)
        totals['second_order'][0] += lo+tail_lo
        totals['second_order'][1] += lo+F(periods, scale)+tail_hi
    totals['total'] = [sum(v[i] for v in totals.values()) for i in (0, 1)]
    assert totals['total'][0] == F(30458547243805493774844991, 10**24)
    return totals


def final_bounds(spectral, density):
    """Conservative rational arithmetic in the proof of Theorem 6.16."""
    g_lower = F('112.8075')
    gamma_lower = F('30.458547')
    raw_upper = F('53.287584')
    combined_upper = F('65.451536565866')
    assert spectral['g'][0] > g_lower
    assert spectral['raw_Q_growth'][1] < raw_upper
    assert spectral['scaled_Q_growth'][1] < combined_upper
    assert density['total'][0] > gamma_lower
    certified_h_upper = combined_upper+57-density['total'][0]
    h_upper = F('91.992990')
    assert certified_h_upper < h_upper
    delta_lower = g_lower/h_upper-1
    mu_upper = g_lower/(g_lower-h_upper)
    assert delta_lower > F('0.2262')
    assert mu_upper < F('5.4197')
    # The archive also certifies a slightly stronger readable delta.
    sharper_h = spectral['scaled_Q_growth'][1]+57-density['total'][0]
    sharper_delta = spectral['g'][0]/sharper_h-1
    assert sharper_delta > F('0.22626234')
    return dict(g_lower=g_lower, h_upper=h_upper, delta_lower=delta_lower,
                mu_upper=mu_upper, sharper_delta_lower=sharper_delta)
