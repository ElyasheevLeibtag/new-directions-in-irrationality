"""Exact prime counts for the zeta(3), zeta(4) proof.

This is a Python port of the supplied conservative prime-count certificate.
It supports Lemma 6.14 (prop:z34-divisor) of the current manuscript.
The function certified_divisor returns the paper's mathcal D_n.
All bounds, geometry, integration, and comparisons use integers/Fractions.
No floating-point value participates in a certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, gcd, isqrt, lcm

INF = 10**9  # Larger than any valuation used here; means 'absent term'.
FORWARD_SLOPES = (-14, -13, -12, -4, -3, -2, -1, 0, 14, 15)
DUAL_SLOPES = (-17, -16, -2, -1, 0, 1, 2, 10, 11, 12)
LAMBDA_EXPONENTS = {8: 1, 9: 2, 11: 3, 12: 3, 13: 2,
                    16: -2, 17: -2, 18: -2, 19: -1}


def decimal_text(q, places=12):
    """Truncated decimal for display, computed by integer long division."""
    q = F(q)
    sign = '-' if q < 0 else ''
    q = abs(q)
    integer, remainder = divmod(q.numerator, q.denominator)
    digits = []
    for _ in range(places):
        digit, remainder = divmod(10*remainder, q.denominator)
        digits.append(str(digit))
    return f'{sign}{integer}.'+''.join(digits)


def valuation(q, p):
    """Exact v_p(q); zero has infinite valuation."""
    q = F(q)
    if not q:
        return float('inf')
    a, b, v = abs(q.numerator), q.denominator, 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def floor_offset(a, d, offset=0, finite=False):
    """Floor(a/d + offset/p): finite d=p, or the p -> infinity limit."""
    if finite:
        return (a + offset) // d
    return a // d - int(a % d == 0 and offset < 0)


def binomial_valuation(a, b, d):
    return a // d - b // d - (a-b) // d


def forward_term_bound(x, y, d, finite=False):
    """x/d=n/p, y/d=k/p. Bound for both numerator readouts."""
    v = lambda a, c=0: floor_offset(a, d, c, finite)
    b = lambda a, c: binomial_valuation(a, c, d)
    h = int(4*x+y >= d)
    common = b(y+12*x, 10*x)+b(y+13*x, 10*x)+b(y+14*x, 10*x)
    if 0 <= y <= 14*x:
        vc = common+b(15*x, x+y)+b(15*x, y)
        slope = vc-1
    elif 14*x < y <= 15*x:
        vc = INF
        slope = common+b(15*x, y)-v(x+y)+v(15*x)+v(y-14*x, -1)
    elif -x <= y < 0:
        vc = INF
        slope = common+b(15*x, x+y)-v(15*x-y)+v(-y, -1)+v(15*x)
    else:
        raise ValueError('forward index outside its support')
    return min(vc-4*h, slope-3*h)


def dual_pole_data(x, y, d, finite=False):
    """Pole order, leading Laurent valuation, and harmonic-layer bound."""
    v = lambda a, c=0: floor_offset(a, d, c, finite)
    shift = int(finite)
    m = int(y <= 11*x)+int(x <= y <= 12*x)+int(2*x <= y <= 10*x+shift)
    t1 = v(y)+v(11*x-y) if y <= 11*x else v(y)-v(y-11*x, -1)
    t2 = v(y-x)+v(12*x-y) if y >= x else v(12*x-y)-v(x-y, -1)
    if 2*x <= y <= 10*x+shift:
        t3 = v(y-2*x)+v(10*x-y, 1)
    elif y < 2*x:
        t3 = v(10*x-y, 1)-v(2*x-y, -1)
    else:
        t3 = v(y-2*x)-v(y-10*x, -2)
    lead = v(16*x+y)-v(2*x+y)+v(17*x+y)-v(x+y)-t1-t2-t3
    harmonic = int(2*x+y >= d)
    layer = min(lead-(m-i)-(i+1)*harmonic for i in range(1, m+1))
    return m, lead, layer


def phi_lambda(x, d):
    return sum(e*(a*x//d) for a, e in LAMBDA_EXPONENTS.items())


def dual_count_samples(x, samples, d, finite=False, phi=None):
    amin = bmin = dmin = INF
    for y in samples:
        m, lead, layer = dual_pole_data(x, y, d, finite)
        if m == 3:
            amin, bmin = min(amin, lead), min(bmin, lead-1)
        elif m == 2:
            bmin = min(bmin, lead)
        dmin = min(dmin, layer)
    if phi is None:
        phi = phi_lambda(x, d)
    return max(0, 4*int(19*x >= d)+phi+min(amin+bmin, amin+dmin, bmin+dmin))


def y_samples(x, d, edges, slopes):
    """All floor breaks, endpoints, and midpoints, represented over 2*d.

    Each structural/harmonic band needs at most one whole y-period.
    The original far endpoints are retained, including their offsets.
    """
    edges = sorted(set(edges))
    samples = {2*y for y in edges}
    for lo, hi in zip(edges, edges[1:]):
        cut = min(hi, lo+d)
        breaks = {lo, cut}
        for slope in slopes:
            first = -((slope*x-lo)//d)  # ceil((lo-slope*x)/d)
            last = (cut-slope*x)//d
            breaks.update(m*d+slope*x for m in range(first, last+1))
        breaks = sorted(breaks)
        samples.update(2*y for y in breaks)
        samples.update(a+b for a, b in zip(breaks, breaks[1:]))
    return samples


@lru_cache(None)
def counts_core(point):
    """Unreduced exact counts at any positive rational x."""
    point = F(point)
    x, d = point.numerator, point.denominator
    if x <= 0:
        raise ValueError('x must be positive')
    edges = [-x, 0, 14*x, 15*x]
    if -x < d-4*x < 15*x:
        edges.append(d-4*x)
    fs = y_samples(x, d, edges, FORWARD_SLOPES)
    forward = max(0, 4*int(19*x >= d)+min(forward_term_bound(2*x, y, 2*d) for y in fs))
    edges = [0, x, 2*x, 10*x, 11*x, 12*x]
    if 0 < d-2*x < 12*x:
        edges.append(d-2*x)
    ds = y_samples(x, d, edges, DUAL_SLOPES)
    dual = dual_count_samples(2*x, ds, 2*d)
    return forward, dual


def counts(point):
    point = F(point)
    if point >= 2:
        point = 1+point-point.numerator//point.denominator
    return counts_core(point)


def rational_grid(order=40):
    return sorted({F(j, d) for d in range(1, order+1) for j in range(d+1)})


def offset_state_counts(point):
    """Include EVERY integer state in thin shifted-wall neighborhoods.

    For x=a/b, y=Y/(2b), set model n=2*a*T, p=2*b*T, k=Y*T+delta.
    These are only integer floor arguments: model p need not be prime.
    T=16 makes every nonzero unshifted distance at least 16, whereas the
    perturbations delta=-3,...,3 and factorial offsets |c|<=2 total at most 5.
    The resulting literal floors therefore equal their exact chamber limits.
    """
    point = F(point)
    a, b = point.numerator, point.denominator
    assert a > 0
    scale = 16
    n, p = 2*a*scale, 2*b*scale
    edges = [-a, 0, 14*a, 15*a]
    if -a < b-4*a < 15*a:
        edges.append(b-4*a)
    samples = y_samples(a, b, edges, FORWARD_SLOPES)
    indices = {y*scale+delta for y in samples for delta in range(-3, 4)
               if -n <= y*scale+delta <= 15*n}
    forward = max(0, 4*int(19*n >= p)
                  +min(forward_term_bound(n, k, p, True) for k in indices))
    edges = [0, a, 2*a, 10*a, 11*a, 12*a]
    if 0 < b-2*a < 12*a:
        edges.append(b-2*a)
    samples = y_samples(a, b, edges, DUAL_SLOPES)
    indices = {y*scale+delta for y in samples for delta in range(-3, 4)
               if 0 <= y*scale+delta <= 12*n}
    dual = dual_count_samples(n, indices, p, True, phi_lambda(n, p))
    return forward, dual


def verify_offset_chambers(order=40):
    """Exact finite chamber certificate, not a sample of integer n values."""
    assert order >= 40
    grid = rational_grid(order)
    checked = 0
    for shift in (0, 1):
        for a, b in zip(grid, grid[1:]):
            point = shift+(a+b)/2
            ordinary = counts_core(point)
            expanded = offset_state_counts(point)
            assert expanded == ordinary, (point, ordinary, expanded)
            if point < F(1, 31):
                assert ordinary == (0, 0)
            checked += 1
    return {'checked_chambers': checked, 'integer_offsets': list(range(-3, 4)),
            'model_scale': 16, 'status': 'PASS'}


def density_certificate(cutoff=40, order=40, diagnostics=True):
    """Exact integrals on (0,cutoff], periodic tail bounds, and cell checks."""
    assert cutoff >= 2 and order >= 40
    grid = rational_grid(order)
    base = []
    for shift in (0, 1):
        rows = []
        for a, b in zip(grid, grid[1:]):
            a, b = a+shift, b+shift
            values = counts_core((a+b)/2)
            if diagnostics:
                assert counts_core((2*a+b)/3) == values
                assert counts_core((a+2*b)/3) == values
            rows.append(values)
        base.append(rows)
    assert base[0][0] == (0, 0)
    integrals = [F(0), F(0), F(0)]
    rows = []
    for shift in range(cutoff):
        for j, (a, b) in enumerate(zip(grid, grid[1:])):
            a, b = a+shift, b+shift
            f, dual = base[int(shift > 0)][j]
            weight = 1/a-1/b if a else F(0)
            for i, value in enumerate((f, dual, max(f, dual))):
                integrals[i] += value*weight
            rows.append((a, b, f, dual, weight))
    assert sum(a*e for a, e in LAMBDA_EXPONENTS.items()) == 0
    assert counts_core(F(1)) == counts_core(F(2))
    for point in map(F, ['7/3', '19/5', '17/2', '159/4', '40']):
        assert counts_core(point) == counts(point)
    endpoint_values = [counts_core(1+x) for x in grid]
    maxima = [max(v[i] for v in base[1]+endpoint_values) for i in (0, 1)]
    tails = [F(m, cutoff) for m in maxima]
    return dict(cutoff=cutoff, cells=len(rows), integrals=integrals,
                maxima=maxima, tail_upper=tails,
                individual_upper=[integrals[i]+tails[i] for i in (0, 1)],
                rows=rows)


def full_density_interval(certificate, tail_start=1000, places=24):
    """Rational interval for the full combined integral.

    Sum the periodic cells exactly with directed fixed-denominator rounding.
    If A is the integral of g_best over one unit period, monotonicity of
    1/x**2 and the integral test bound the remaining tail beyond integer J by
        A/(J+1) <= tail <= A*(1/J + 1/J**2).
    No digamma evaluation or floating-point error estimate is used.
    """
    start = certificate['cutoff']
    assert tail_start > start >= 2 and places >= 12
    period = [(a-1, b-1, max(f, d)) for a, b, f, d, _ in certificate['rows']
              if 1 <= a and b <= 2]
    mean = sum((v*(b-a) for a, b, v in period), F(0))
    scale = 10**places
    lower_units = upper_units = 0
    for j in range(start, tail_start):
        for a, b, value in period:
            if value == 0:
                continue
            # value*(1/(j+a)-1/(j+b)) as an exact numerator/denominator.
            numerator = value*(b.numerator*a.denominator-a.numerator*b.denominator)
            denominator = (j*a.denominator+a.numerator)*(j*b.denominator+b.numerator)
            quotient, remainder = divmod(scale*numerator, denominator)
            lower_units += quotient
            upper_units += quotient+int(remainder != 0)
    lower_tail = mean/F(tail_start+1)
    upper_tail = mean*(F(1, tail_start)+F(1, tail_start**2))
    lower = certificate['integrals'][2]+F(lower_units, scale)+lower_tail
    upper = certificate['integrals'][2]+F(upper_units, scale)+upper_tail
    return dict(lower=lower, upper=upper, period_mean=mean,
                tail_start=tail_start, directed_rounding_places=places,
                residual_lower=lower_tail, residual_upper=upper_tail)


@lru_cache(None)
def harmonic(n, order=1):
    return sum((F(1, j**order) for j in range(1, n+1)), F(0))


def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(None)
def forward_B(n, k):
    return (binom(k+12*n, 10*n)*binom(k+13*n, 10*n)*binom(k+14*n, 10*n)
            *binom(15*n, n+k)*binom(15*n, k))


@lru_cache(None)
def forward_slope(n, k):
    if 0 <= k <= 14*n:
        slope = sum(harmonic(k+(10+u)*n)-harmonic(k+u*n) for u in (2, 3, 4))
        slope += harmonic(14*n-k)-harmonic(n+k)+harmonic(15*n-k)-harmonic(k)
        return forward_B(n, k)*slope
    if 14*n < k <= 15*n:
        num = (binom(k+12*n, 10*n)*binom(k+13*n, 10*n)*binom(k+14*n, 10*n)
               *binom(15*n, k))
        return F((-1)**k*num, (n+k)*binom(n+k-1, 15*n))
    if -n <= k < 0:
        i = -k
        num = (binom(12*n-i, 10*n)*binom(13*n-i, 10*n)*binom(14*n-i, 10*n)
               *binom(15*n, n-i))
        return F((-1)**(i-1)*num, (15*n+i)*binom(15*n+i-1, i-1))
    return F(0)


def exact_readouts(n):
    q = sum(forward_B(n, k) for k in range(14*n+1))
    numerators = [sum((forward_B(n, k)*harmonic(4*n+k, s)
                      -forward_slope(n, k)*harmonic(4*n+k, s-1)/(s-1)
                      for k in range(-n, 15*n+1)), F(0)) for s in (3, 4)]
    return [F(q), *numerators]


def lambda_scalar(n):
    result = F(1, 6)
    for a, exponent in LAMBDA_EXPONENTS.items():
        result *= F(factorial(a*n))**exponent
    return result


def primes_below(limit):
    return [p for p in range(2, limit) if all(p % q for q in range(2, isqrt(p)+1))]


def certified_divisor(n):
    """Construct the paper's common divisor mathcal D_n from floor bounds.

    The returned tuple is (divisor, {prime: exponent}). It never computes
    the actual gcd; the notebook compares that gcd only as a diagnostic.
    """
    assert isinstance(n, int) and n >= 1
    lam = lambda_scalar(n)
    factors = {}
    divisor = 1
    for p in primes_below(31*n):
        if p <= 3 or p*p <= 30*n:
            continue
        forward = max(0, 4*int(p <= 19*n)
                      +min(forward_term_bound(n, k, p, True) for k in range(-n, 15*n+1)))
        dual = dual_count_samples(n, range(12*n+1), p, True, valuation(lam, p))
        exponent = max(forward, dual)
        if exponent:
            factors[p] = exponent
            divisor *= p**exponent
    return divisor, factors


def finite_prime_checks(indices=range(4, 9)):
    """Literal integer floors, including offsets. Independent exact gcd checks."""
    rows, captures = [], {}
    for n in indices:
        multiplier = 6*lcm(*range(1, 19*n+1))**4
        scaled = [multiplier*q for q in exact_readouts(n)]
        assert all(q.denominator == 1 for q in scaled)
        actual_gcd = gcd(*(q.numerator for q in scaled))
        lam = lambda_scalar(n)
        captures[n] = []
        for p in primes_below(31*n):
            if p*p <= 30*n:
                continue
            forward = max(0, 4*int(p <= 19*n)
                          +min(forward_term_bound(n, k, p, True) for k in range(-n, 15*n+1)))
            dual = dual_count_samples(n, range(12*n+1), p, True, valuation(lam, p))
            actual = valuation(actual_gcd, p)
            assert max(forward, dual) <= actual, (n, p, forward, dual, actual)
            rows.append((n, p, forward, dual, actual))
            if dual == actual and dual > forward:
                captures[n].append(p)
    return rows, captures
