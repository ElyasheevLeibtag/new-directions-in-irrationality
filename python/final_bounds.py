"""Exact rational arithmetic for the paper's three reported consequences.

This module combines rate bounds; it does not prove their analytic or
arithmetic hypotheses. Those computations and the remaining mathematical
inputs are described in the subject notebooks.
"""
from fractions import Fraction

# g is a lower bound for the exponential approximation-error rate.
# h is an upper bound for the growth of an integer denominator representative.
# The reduced denominator is no larger than this representative.
PAPER_BOUNDS = {
    'A*zeta(3)+B*zeta(4)': {
        'g': '112.8075', 'h': '111.105', 'delta_claim': '0.015',
        'mu_claim': '66.3', 'label': 'thm:intro-zeta34',
    },
    'zeta(2)': {
        'g': '62.636575', 'h': '49.098639', 'delta_claim': '0.2757',
        'mu_claim': '4.62675', 'label': 'thm:z23-measures',
    },
    'zeta(3)': {
        'g': '112.8075', 'h': '91.992990', 'delta_claim': '0.2262',
        'mu_claim': '5.4197', 'label': 'thm:z23-measures',
    },
}

def consequences(g, h):
    """Return g/h-1 and g/(g-h), using decimal strings or exact rationals."""
    g, h = Fraction(g), Fraction(h)
    if not (g > h > 0):
        raise ValueError('Positive approximation exponent requires g > h > 0.')
    return g/h-1, g/(g-h)

def verify_final_rounding():
    rows = []
    for target, data in PAPER_BOUNDS.items():
        delta, mu = consequences(data['g'], data['h'])
        assert delta > Fraction(data['delta_claim'])
        assert mu < Fraction(data['mu_claim'])
        rows.append({'target': target, **data, 'delta_lower': str(delta),
                     'mu_upper': str(mu), 'delta_decimal': float(delta),
                     'mu_decimal': float(mu)})
    return rows

if __name__ == '__main__':
    import json
    print(json.dumps(verify_final_rounding(), indent=2))
