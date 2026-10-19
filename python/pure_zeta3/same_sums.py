"""Independent finite checks of the all-index identity in Appendix B.4.

The proof is the gamma-reflection/Laurent-coefficient argument in notebook 01.
Here the forward binomial interpolation and its two simple-zero tails are
computed independently of the native kernel partial-fraction implementation.
"""
from math import comb, factorial
from flint import fmpq
from approximants import partial_fractions_direct, kernel_readout


def harmonic(m, order=1):
    return sum((fmpq(1, j**order) for j in range(1, m+1)), fmpq(0))


def forward_coefficients(n, k):
    """Return B_{n,k}, B'_{n,k}, including -n <= k <= 15n tails."""
    if n < 0 or not -n <= k <= 15*n:
        raise ValueError('Require n >= 0 and -n <= k <= 15n.')
    polynomial = comb(k+12*n, 10*n)*comb(k+13*n, 10*n)*comb(k+14*n, 10*n)
    if 0 <= k <= 14*n:
        B = fmpq(polynomial*comb(15*n, n+k)*comb(15*n, k))
        derivative = sum((harmonic(k+(10+u)*n)-harmonic(k+u*n)
                          for u in (2, 3, 4)), fmpq(0))
        derivative += harmonic(14*n-k)-harmonic(n+k)+harmonic(15*n-k)-harmonic(k)
        return B, B*derivative
    if k < 0:
        i = -k
        derivative = fmpq((-1)**(i-1)*polynomial*comb(15*n,n-i),
                          (15*n+i)*comb(15*n+i-1,i-1))
    else:
        derivative = fmpq((-1)**k*polynomial*comb(15*n,k),
                          (n+k)*comb(n+k-1,15*n))
    return fmpq(0), derivative


def check_same_sums(n):
    """Check every pole and both full sums at a specified finite index."""
    Cprime = fmpq(factorial(15*n)**2, factorial(10*n)**3)
    coefficients = partial_fractions_direct(n)
    bold_q, bold_p3 = fmpq(0), fmpq(0)
    for k in range(-n, 15*n+1):
        B, Bprime = forward_coefficients(n,k)
        a, b = coefficients[2*n+k]
        assert Cprime*a == (-1)**n*B
        assert Cprime*b == -(-1)**n*Bprime
        bold_q += B
        bold_p3 += B*harmonic(4*n+k,3)-Bprime*harmonic(4*n+k,2)/2
    Q, P3, _ = kernel_readout(n)
    assert Cprime*Q == (-1)**n*bold_q
    assert Cprime*P3 == (-1)**n*bold_p3
    assert bold_q > 0 and P3/Q == bold_p3/bold_q
    return {'n':n, 'coefficient_pairs':len(coefficients), 'same_ratio':True}
