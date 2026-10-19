"""Exact reconstruction of Lemma 6.17 (lem:z2-connection).

The full M(n), Q(n), K(n), c(n) are derived from the supplied manuscript's
parameter data, using RamanujanTools. No unavailable archive file is required.
"""
import sympy as sp
from ramanujantools import Matrix, Position
from ramanujantools.cmf import pFq, MeijerG
from ramanujantools.flint_core import SymbolicMatrix, flint_ctx


def certify_pairing():
    """Return the rational pairing after verifying M K(n+1) Q^T = c K.

    C is the bilinear concomitant at z=1 for P(θ)-zR(θ). The direct
    adjoint Meijer function has upper parameters α. The paper increases
    the last upper parameter by one. Its contiguous relation gives
    H_state = G_state A, A=α_3 I-Θ_G. Therefore K=C A^T.

    The concomitant derivative formula below includes derivatives of the
    coefficient z: at z=1 every positive θ-derivative of -z R_i is -R_i,
    while a zeroth derivative is P_i-R_i. The formula is computed in the
    derivative basis (F,θF,θ²F).
    """
    n, theta = sp.symbols("n theta")
    context = flint_ctx({n}, fmpz=True)

    def rational_matrix(matrix):
        return SymbolicMatrix.from_sympy(Matrix(matrix), context)

    alpha = [-11*n, 1+10*n, 1+9*n, 1+11*n]
    beta = [1+2*n, 1+n, 1+3*n]
    P = sp.Poly(theta*sp.prod(theta+y-1 for y in beta), theta)
    R = sp.Poly(sp.prod(theta+x for x in alpha), theta)
    concomitant = sp.zeros(3)
    for row in range(3):
        for order in range(row+1, 5):
            derivative_order = order-1-row
            for column in range(min(derivative_order, 2)+1):
                coefficient = (P.nth(order)-R.nth(order)
                               if derivative_order == column else -R.nth(order))
                concomitant[row, column] += (
                    (-1)**derivative_order
                    * sp.binomial(derivative_order, column)*coefficient)

    meijer_upper = alpha[:-1]+[alpha[-1]+1]
    meijer_lower = [0]+[y-1 for y in beta]
    meijer_equation = sp.Poly(
        sp.prod(theta-a+1 for a in meijer_upper)
        - sp.prod(theta-b for b in meijer_lower), theta)
    theta_g = Matrix.companion(sp.Poly(
        meijer_equation.as_expr()/meijer_equation.LC(), theta))
    contiguous_shift = alpha[-1]*sp.eye(3)-theta_g
    K = Matrix((concomitant*contiguous_shift.T).applyfunc(sp.factor))

    forward_axes = sp.symbols("x0:4")+sp.symbols("y0:3")
    forward_start = Position(dict(zip(forward_axes, [0,1,1,1,1,1,1])))
    forward_direction = Position(dict(zip(forward_axes, [-11,10,9,11,2,1,3])))
    M = pFq(4,3,1).trajectory_matrix(forward_direction, forward_start, n)

    dual_axes = sp.symbols("a0:4")+sp.symbols("b0:4")
    dual_start = Position(dict(zip(dual_axes, [0,1,1,2,0,0,0,0])))
    dual_direction = Position(dict(zip(dual_axes, [-11,10,9,11,0,2,1,3])))
    Q = MeijerG(4,1,4,4,-1).trajectory_matrix(dual_direction, dual_start, n)

    product = (rational_matrix(M)*rational_matrix(K).shift(n,1)
               * rational_matrix(Q.T)).factor()
    scalar = sp.factor(product[0,2]/K[0,2])
    assert rational_matrix(product) == rational_matrix(scalar*K)
    expected_determinant = -sp.prod(j*n+1 for j in (8,9,10,11,13))*(13*n+2)
    assert sp.factor(K.det()-expected_determinant) == 0
    assert K.subs({n:0}) == Matrix([[0,0,-1],[-1,2,-1],[-2,2,0]])
    return {"n":n,"M":M,"Q":Q,"K":K,"c":scalar,
            "determinant_K":sp.factor(K.det()),"identity_verified":True}
