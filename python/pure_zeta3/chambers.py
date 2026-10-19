"""Exact arithmetic certificates for Appendix ``app:zeta3-measure``.

The algorithms originate in the supplied archive's exact_chambers.py,
all_shell_cells.py, high_kernel.py, and second_order_certificate.py.
Only exact Fraction arithmetic is used. A midpoint represents an entire open
cell because the listed affine floor boundaries cannot cross inside it;
all y-boundary strata are also checked. The notebook explains the distinct
uniform analytic/divisibility inputs needed to interpret these certificates.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
from math import floor, ceil
import json
from pathlib import Path

def floor_fraction(x):
    """Exact floor; unlike int(), this is correct for negative rational inputs."""
    return x.numerator // x.denominator

def ceil_fraction(x):
    """Exact ceiling of a rational number."""
    return -floor_fraction(-x)

@lru_cache(None)
def harmonic_number(n):
    """The exact ordinary harmonic number H_n, with H_0=0."""
    return sum((F(1, k) for k in range(1, n + 1)), F(0))

def coefficient_valuation_bounds(x, y, kind):
    """Termwise lower bounds for the order-3 and order-4 normalized readouts.

x=n/p, y=j/p. kind='d' is a double pole, 'l' a left simple
pole, and 'r' a right simple pole. The exact harmonic coefficient
in the logarithmic derivative decides whether a prime is lost."""
    c = 2 * floor_fraction(15 * x) - 3 * floor_fraction(10 * x)
    num = c + sum((floor_fraction(y + (10 + u) * x) - floor_fraction(y + u * x) for u in (0, 1, 2)))
    h = int(y + 2 * x >= 1)
    l = int(x >= F(1, 19))
    if kind == 'd':
        s = num - floor_fraction(y - x) - floor_fraction(y - 2 * x) - floor_fraction(16 * x - y) - floor_fraction(17 * x - y)
        lp = (
            sum(
                -harmonic_number(floor_fraction(y + (10 + u) * x))
                + harmonic_number(floor_fraction(y + u * x))
                for u in (0, 1, 2)
            )
            + harmonic_number(floor_fraction(y - x))
            - harmonic_number(floor_fraction(16 * x - y))
            + harmonic_number(floor_fraction(y - 2 * x))
            - harmonic_number(floor_fraction(17 * x - y))
        )
        loss = int(lp != 0)
        return [min(s + r * l - r * h, s - loss + r * l - (r - 1) * h) for r in (3, 4)]
    if kind == 'l':
        b = num + ceil_fraction(2 * x - y) - 1 - floor_fraction(16 * x - y) - floor_fraction(y - x) - floor_fraction(17 * x - y)
    else:
        b = num + ceil_fraction(y - 16 * x) - 1 - floor_fraction(y - 2 * x) - floor_fraction(y - x) - floor_fraction(17 * x - y)
    return [b + r * l - (r - 1) * h for r in (3, 4)]

def baseline_profile(x):
    """Minimize the readout valuations on every y-cell and boundary.

The affine breakpoints cover every factorial floor and harmonic
cutoff. This computes g3 and g4 in the first period, or their
common periodic tail when x is shifted by 1."""
    gains = [100, 100]
    for (kind, lo, hi, include_lo, include_hi) in [('d', 2 * x, 16 * x, True, True), ('l', x, 2 * x, True, False), ('r', 16 * x, 17 * x, False, True)]:
        cuts = {lo, hi}
        for slope in (-12, -11, -10, -2, -1, 0, 1, 2, 16, 17):
            for k in range(floor_fraction(lo - slope * x) - 1, ceil_fraction(hi - slope * x) + 2):
                v = k + slope * x
                if lo <= v <= hi:
                    cuts.add(v)
        cuts = sorted(cuts)
        samples = [(a + b) / 2 for (a, b) in zip(cuts, cuts[1:])]
        samples += [y for y in cuts if (y != lo or include_lo) and (y != hi or include_hi)]
        for y in samples:
            gs = coefficient_valuation_bounds(x, y, kind)
            gains = [min(a, b) for (a, b) in zip(gains, gs)]
    return [max(0, g) for g in gains]

@lru_cache(None)
def high_kernel_readouts(numerator_roots, denominator_roots, harmonic_shift):
    """Evaluate the local partial-fraction readouts using exact rationals.

numerator_roots/denominator_roots are integer high-root multisets;
harmonic_shift is the shifted cutoff. Equal output on all minimal
classes is the first-order shell criterion."""
    H = lambda m, r: sum((F(1, k ** r) for k in range(1, m + 1)), F(0))
    ks = []
    residues = []
    for (l, m) in Counter(denominator_roots).items():
        f = F(1)
        for z in numerator_roots:
            f *= l - z
        for z in denominator_roots:
            if z != l:
                f /= l - z
        if m == 1:
            a = F(0)
            b = f
        elif m == 2:
            a = f
            b = f * (sum((F(1, l - z) for z in numerator_roots), F(0)) - sum((F(1, l - z) for z in denominator_roots if z != l), F(0)))
        else:
            raise ValueError((numerator_roots, denominator_roots, harmonic_shift))
        residues.append((l, a, b))
    for r in (3, 4):
        z = sum((a * H(-l + harmonic_shift, r) + b * H(-l + harmonic_shift, r - 1) / (r - 1) for (l, a, b) in residues), F(0))
        ks.append(str(z))
    return ks
SLOPES = (0, 1, 2, 10, 11, 12, -1, -2, -16, -17)

def high_root_signature(x, y):
    """Return high roots and cutoff for x=n/p and y=k/p.

The explicit floor/ceiling intervals include boundary residues.
Only open x-cells are used because eligible primes avoid Farey walls."""
    numerator_roots = []
    denominator_roots = []
    for u in (0, 1, 2):
        numerator_roots += [l for l in range(floor(u * x - y) + 1, floor((10 + u) * x - y) + 1)]
    for (a, b) in ((1, 17), (2, 16)):
        denominator_roots += [l for l in range(ceil(-b * x - y), floor(-a * x - y) + 1)]
    return (tuple(sorted(numerator_roots)), tuple(sorted(denominator_roots)), floor(2 * x - y))

def first_order_cell(a, b):
    """Verify the complete first-order root/kernel criterion on one x-cell.

Each y-interval midpoint and each y-boundary is visited. The returned
K values are exact rationals, not sampled numerical estimates."""
    x = (a + b) / 2
    cuts = {F(0), F(1)}
    for m in SLOPES:
        for c in range(-20, 21):
            y = m * x + c
            if 0 <= y <= 1:
                cuts.add(y)
    cuts = sorted(cuts)
    points = cuts[:-1] + [(u + v) / 2 for (u, v) in zip(cuts, cuts[1:])]
    signatures = set()
    normalized_endpoints = set()
    hmin = 0
    for y in points:
        (numerator_roots, denominator_roots, harmonic_shift) = high_root_signature(x, y)
        h = len(numerator_roots) - len(denominator_roots)
        hmin = min(hmin, h)
        if h == -1:
            signatures.add((numerator_roots, denominator_roots, harmonic_shift))
            z = []
            for u in (0, 1, 2):
                z.extend([floor(u * x - y) + 1 - harmonic_shift, floor((u + 10) * x - y) - harmonic_shift])
            for (aa, bb) in ((1, 17), (2, 16)):
                z.extend([ceil(-bb * x - y) - harmonic_shift, floor(-aa * x - y) - harmonic_shift])
            normalized_endpoints.add(tuple(z))
    if hmin < -1:
        return dict(valid=False, hmin=hmin, point_count=len(points))
    ks = {tuple(high_kernel_readouts(*z)) for z in signatures}
    if len(ks) != 1:
        return dict(valid=False, hmin=hmin, K_count=len(ks), point_count=len(points))
    assert len(normalized_endpoints) == 1, (a, b, normalized_endpoints)
    return dict(
        valid=True,
        normalized_interval_endpoints=list(next(iter(normalized_endpoints))),
        hmin=hmin,
        K=list(next(iter(ks))),
        point_count=len(points),
        signature_count=len(signatures),
        signatures=[
            dict(ns=numerator_roots, ds=denominator_roots, eps=harmonic_shift)
            for numerator_roots, denominator_roots, harmonic_shift in sorted(signatures)
        ],
        C_v=2 * floor(15 * x) - 3 * floor(10 * x),
    )

def affine_root_intervals(x, y):
    """Store high-root intervals with affine endpoints (slope, intercept).

An endpoint (m,c) means m*t+c for the symbolic periodic lift t.
This representation proves all lifts at once."""
    numerator_roots = []
    denominator_roots = []
    for u in (0, 1, 2):
        numerator_roots.append(((u, floor(u * x - y) + 1), (u + 10, floor((u + 10) * x - y))))
    for (a, b) in ((1, 17), (2, 16)):
        denominator_roots.append(((-b, ceil(-b * x - y)), (-a, floor(-a * x - y))))
    return (tuple(numerator_roots), tuple(denominator_roots))

def affine_quotient_factor(cur, ref):
    """Certify that current/reference is exactly one factor V-a(t).

Signed interval boundaries cancel algebraically; two remaining
opposite boundaries at a(t) and a(t)+1 encode a single root."""
    c = Counter()
    for (sign, ints) in ((1, cur[0]), (1, ref[1]), (-1, cur[1]), (-1, ref[0])):
        for (a, b) in ints:
            c[a] += sign
            c[b[0], b[1] + 1] -= sign
    c = {k: v for (k, v) in c.items() if v}
    pos = [k for (k, v) in c.items() if v == 1]
    neg = [k for (k, v) in c.items() if v == -1]
    if len(c) == 2 and len(pos) == len(neg) == 1 and (neg[0] == (pos[0][0], pos[0][1] + 1)):
        return pos[0]
    return None

def second_order_cell(a, b):
    """Check all affine trace-lemma hypotheses on one Farey x-cell.

Minimal classes must have one common high-root function. Every
zero class must differ by one affine factor, and the residue-zero
multiplicity must be positive."""
    x = (a + b) / 2
    cuts = {F(0), F(1)}
    for m in SLOPES:
        for c in range(-20, 21):
            y = m * x + c
            if 0 <= y <= 1:
                cuts.add(y)
    cuts = sorted(cuts)
    points = cuts[:-1] + [(u + v) / 2 for (u, v) in zip(cuts, cuts[1:])]
    minimal_signatures = set()
    zeros = set()
    hmin = 0
    for y in points:
        (numerator_roots, denominator_roots, e) = high_root_signature(x, y)
        h = len(numerator_roots) - len(denominator_roots)
        hmin = min(hmin, h)
        if h == -1:
            minimal_signatures.add(affine_root_intervals(x, y))
        if h == 0:
            zeros.add(affine_root_intervals(x, y))
    if hmin < -1 or len(minimal_signatures) != 1:
        return dict(valid=False, hmin=hmin, min_signatures=len(minimal_signatures))
    ref = next(iter(minimal_signatures))
    quotients = [affine_quotient_factor(z, ref) for z in zeros]
    if any((q is None for q in quotients)):
        return dict(valid=False, hmin=hmin, min_signatures=1, failed_h0=True)
    (numerator_roots, denominator_roots, e) = high_root_signature(x, F(0))
    h0 = len(numerator_roots) - len(denominator_roots)
    if h0 < 1:
        return dict(valid=False, hmin=hmin, residue_zero_h=h0)
    return dict(valid=True, hmin=hmin, min_signature=ref, h0_affine_factors=sorted(set(quotients)), residue_zero_h=h0, point_count=len(points))

def build_chambers():
    """Recompute all 270 Farey-29 cells and both shell criteria."""
    endpoints = sorted({F(a, b) for b in range(1, 30) for a in range(b + 1)})
    rows = []
    for (left, right) in zip(endpoints, endpoints[1:]):
        x = (left + right) / 2
        (first_period, periodic_tail) = (baseline_profile(x), baseline_profile(x + 1))
        assert periodic_tail[0] == periodic_tail[1]
        row = dict(lo=str(left), hi=str(right), g3=first_period[0], g4=first_period[1], tail=periodic_tail[0])
        row.update(first_order_cell(left, right))
        if row['valid']:
            assert row['C_v'] + 2 - row['tail'] == 1
            assert row['C_v'] + 2 - row['g3'] == 1
            assert row['C_v'] + 2 - row['g4'] == 1
        row['second'] = second_order_cell(left, right)
        if row['second']['valid']:
            assert row['valid']
            assert row['second']['hmin'] == -1
            assert row['second']['residue_zero_h'] >= 1
        rows.append(row)
    assert len(rows) == 270
    return rows

def merge_runs(rows, second=False):
    """Union adjacent valid cells into the open shell/window intervals."""
    runs = []
    for row in rows:
        valid = row['second']['valid'] if second else row['valid']
        if valid:
            (a, b) = (F(row['lo']), F(row['hi']))
            if runs and runs[-1][1] == a:
                runs[-1] = (runs[-1][0], b)
            else:
                runs.append((a, b))
    return runs

def compare_reference(rows):
    """Compare regenerated rational data with the supplied exact fixtures."""
    reference = json.loads((Path(__file__).parent / 'data/chamber_reference.json').read_text())
    for (actual, base, first, second) in zip(rows, reference['baseline'], reference['first_order'], reference['second_order']):
        for key in ('lo', 'hi', 'g3', 'g4', 'tail'):
            assert actual[key] == base[key], (key, actual)
        for key in ('valid', 'hmin', 'point_count'):
            assert actual[key] == first[key], (key, actual)
        if actual['valid']:
            for key in ('normalized_interval_endpoints', 'K', 'C_v'):
                assert actual[key] == first[key], (key, actual)
        assert json.loads(json.dumps(actual['second'])) == second['second']
    return True
