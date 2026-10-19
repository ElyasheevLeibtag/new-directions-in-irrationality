"""Verify the 5F4 forward limit through the Meijer G dual limit.

Run from the repository root:
    python additional/analytic/verify_5f4_meijerg_recessive.py
"""

import sys
from itertools import permutations

from sympy import Poly, oo, prod, symbols
from sympy.abc import n

from ramanujantools import Matrix, Position
from ramanujantools.cmf import MeijerG, pFq
from ramanujantools.flint_core.context import flint_ctx, flint_to_sympy
from ramanujantools.flint_core.rational import FlintRational
from ramanujantools.flint_core.symbolic_matrix import SymbolicMatrix


def construct_trajectories():
    axesF = symbols("x0:5") + symbols("y0:4")
    startF = Position(dict(zip(axesF, (0, 0, 1, 1, 1, 1, 1, 1, 1), strict=True)))
    directionF = Position(dict(zip(axesF, (-14, -15, 13, 14, 12, 1, 2, 3, 4), strict=True)))
    F = pFq(5, 4, 1).trajectory_matrix(directionF, startF)

    axesG = symbols("a0:5") + symbols("b0:5")
    startG = Position(dict(zip(axesG, (0, 0, 1, 1, 2, 0, 0, 0, 0, 0), strict=True)))
    directionG = Position(dict(zip(axesG, (-16, -17, 11, 12, 10, -2, -1, 0, 1, 2), strict=True)))
    G = MeijerG(5, 2, 5, 5, 1).trajectory_matrix(directionG, startG)

    # Invariant basis on the forward side; L is its constant left inverse.
    E = Matrix([
        [1, 0, 0],
        [0, -n*(450*n + 29), -n**2*(101250*n**2 + 12975*n + 406)],
        [0, 1, 0],
        [0, 0, 1],
    ])
    L = Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    H = Matrix([
        [-6*n**2*(10*n + 1)*(12*n + 1)*(193*n + 13),
         -n*(12*n + 1)*(2550*n**2 + 523*n + 23),
         (12*n + 1)*(450*n**2 + 41*n + 1)],
        [(10*n + 1)*(450*n**2 + 21*n - 1), (-25*n - 2)*(18*n - 1), -1],
        [-20*n - 2, 2, 0],
    ])
    numerator = (
        1125*(n + 1)**3*(2*n + 1)*(3*n + 1)*(3*n + 2)
        * prod(5*n + j for j in range(1, 5))
        * prod(15*n + j for j in (1, 2, 4, 7, 8, 11, 13, 14))
    )
    denominator = (
        416*(6*n + 1)*(6*n + 5)
        * prod(12*n + j for j in (1, 5, 7, 11))
        * prod(13*n + j for j in range(1, 13))
    )
    R = numerator / denominator
    return F, G, E, L, H, R


def determinant3(matrix):
    """Compute an exact 3x3 determinant without a large SymPy expansion."""
    matrix._normalize()
    context = matrix.ctx
    numerator = context.constant(0)
    for permutation in permutations(range(3)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(3) for j in range(i + 1, 3))
        term = context.constant((-1)**inversions)
        for i, j in enumerate(permutation):
            term *= matrix.polynomials[3*i + j]
        numerator += term
    value = FlintRational(numerator, context.constant(1), context)
    for _ in range(3):
        value *= matrix.scalar
    return value


def certify_no_nonnegative_zeros(polynomial):
    """A factor with one coefficient sign and nonzero constant has no n>=0 zero."""
    unit, factors = polynomial.factor()
    assert unit != 0
    for factor, _ in factors:
        p = Poly(flint_to_sympy(factor, factor=False), n)
        assert p.nth(0) != 0
        coefficients = p.all_coeffs()
        assert (all(c >= 0 for c in coefficients)
                or all(c <= 0 for c in coefficients))


def main():
    sys.set_int_max_str_digits(0)
    F, G, E, L, H, R = construct_trajectories()
    context = flint_ctx([n], fmpz=True)
    convert = lambda matrix: SymbolicMatrix.from_sympy(matrix, context)
    Q = G[:3, :3]
    F_restricted = (convert(L)*convert(F)*convert(E.subs({n: n + 1}))).to_rt(factor=True)
    assert L*E == Matrix.eye(3)
    assert G[3, :].factor() == Matrix([[0, 0, 0, -1]])
    assert convert(F)*convert(E.subs({n: n + 1})) == convert(E)*convert(F_restricted)
    assert convert(F_restricted)*convert(H.subs({n: n + 1}))*convert(Q.T) == convert(R*H)
    rectangular_K = E*H
    assert rectangular_K.shape == (4, 3) and rectangular_K.rank() == 3
    assert (convert(F)*convert(rectangular_K.subs({n: n + 1}))*convert(Q.T)
            == convert(R*rectangular_K))
    D = prod(j*n + 1 for j in range(8, 13))
    assert H.det().factor() == 2*D
    assert H.inverse()*Matrix([1, 0, 0]) == Matrix([1, 10*n + 1, (10*n + 1)**2])/D
    print("Invariant restriction and exact coboundary:", True)
    print("Rectangular 4x3 intertwining identity:", True)

    # This certifies that three consecutive error numerators cannot all
    # vanish. Eventual positivity of every mixed error also needs the
    # common-denominator estimate stated in the manuscript.
    reduced = convert(F_restricted)
    reduced_next = convert(F_restricted.subs({n: n + 1}))
    e0 = Matrix([1, 0, 0])
    e0_exact = convert(e0)
    consecutive = convert(Matrix.hstack(
        e0,
        (reduced*e0_exact).to_rt(factor=True),
        (reduced*reduced_next*e0_exact).to_rt(factor=True),
    ))
    for matrix in (reduced, consecutive):
        value = determinant3(matrix)
        certify_no_nonnegative_zeros(value.numerator)
        certify_no_nonnegative_zeros(value.denominator)
    coefficient_map = Matrix([[12, 0, 0], [0, -6, 3], [0, 0, -2]])
    assert coefficient_map.det() == 144
    print("Three consecutive coefficient vectors are independent for every n>=0:", True)

    # K is an intertwining matrix of rank three, not an invertible
    # coboundary of the full four-dimensional systems. The inverse-free
    # identity below implies the exact walk identity for every N.
    J = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    K = E*H*J
    assert K.rank() == 3
    assert J*G.T == Q.T*J
    assert convert(F)*convert(K.subs({n: n + 1}))*convert(G.T) == convert(R*K)
    terminal = Matrix([1, 10*n + 1, (10*n + 1)**2, 0])/D
    assert K*terminal == Matrix([1, 0, 0, 0])
    # Check the telescoped identity with exact finite walks as well.
    WF, WG, scalar = Matrix.eye(4), Matrix.eye(4), 1
    for step in range(3):
        WF = WF*F.subs({n: step})
        WG = WG*G.subs({n: step})
        scalar *= R.subs(n, step)
        assert WF*K.subs({n: step + 1}) == scalar*K.subs({n: 0})*WG.inverse().T
    print("Rank-three intertwining and exact full walk identity:", True)

    A, B, Z2, Z3, Z4 = symbols("A B Z2 Z3 Z4")
    stateG0 = Matrix([[1 - 2*Z3 + 3*Z4, 1 - 2*Z3, 1, 1 - Z2]])
    initialG = Matrix([[2*B, -3*A - 2*B, 3*A, 0], [0, 0, 6, 0]])
    initialF = Matrix([[0, 0, -6*A, 3*A - 2*B], [12, 0, 0, 0]])
    assert initialG*stateG0.T == 6*Matrix([A*Z3 + B*Z4, 1])
    assert initialF*E.subs({n: 0})*H.subs({n: 0}) == 2*initialG[:, :3]
    assert initialF*K.subs({n: 0}) == 2*initialG
    assert initialG*J.T*J == initialG
    # Rows with fourth coordinate zero stay in that space under G^-T:
    # this follows equivalently from J*G.T == Q.T*J above.
    projection = Matrix([1, 0, 0, 0])
    assert E*Matrix([1, 0, 0]) == projection
    assert initialF*K.subs({n: 0})*stateG0.T == 12*Matrix([A*Z3 + B*Z4, 1])
    print("Meijer G initial identity and transport to 5F4:", True)
    print("Initial conditions exclude the -1 mode; denominator pairing is 12:", True)

    # At infinity the transported projection has direction (1,10,100)
    # after balancing the three state coordinates by (1,n^-1,n^-2).
    Q_infinity = Matrix(3, 3, lambda i, j: (n**i*Q[i, j]/(n + 1)**j).limit(n, oo))
    dual_infinity = Q_infinity.inverse().T
    terminal_direction = Matrix([1, 10, 100])
    # The scalar D cancels from the ratio. Balance the terminal column
    # by diag(1,n^-1,n^-2) to obtain its actual limiting direction.
    assert Matrix([(D*terminal[i]/n**i).limit(n, oo) for i in range(3)]) == terminal_direction
    assert Q_infinity.charpoly().discriminant() != 0
    assert Matrix.hstack(terminal_direction, dual_infinity*terminal_direction,
                         dual_infinity**2*terminal_direction).det() != 0
    print("Projection vector couples to every retained limiting mode:", True)

    polynomial = G.charpoly(poincare=True)
    reduced_polynomial = Q.charpoly(poincare=True)
    G_infinity = Matrix(4, 4, lambda i, j: (n**i*G[i, j]/(n + 1)**j).limit(n, oo))
    assert all(entry.is_Rational for entry in G_infinity)
    assert G_infinity.charpoly().all_coeffs() == polynomial.all_coeffs()
    assert polynomial.nth(0) != 0 and polynomial.discriminant() != 0
    assert polynomial.as_expr() == (
        (polynomial.gen + 1)*reduced_polynomial.as_expr()
    ).expand()
    rhos = sorted(
        (root.evalf(40) for root in polynomial.all_roots(radicals=False)),
        key=lambda root: abs(complex(root)),
        reverse=True,
    )
    for i, rho in enumerate(rhos, 1):
        print(f"rho_{i} = {rho}")
    print("Exponential error factor:", (abs(rhos[-1])/abs(rhos[0])).evalf(4))


if __name__ == "__main__":
    main()
