"""Regenerate algebraic and threshold checks for Appendix C.

The printed trajectory is transcribed explicitly. These checks verify its
determinant, limit, characteristic polynomial, cyclic determinant, root-bound
polynomials, and the positive margin at M=7165. They do not replace the
analytical or all-index arithmetic proofs in the manuscript.
"""
from pathlib import Path
import json
import sympy as s
from flint import arb, ctx

if not __debug__:
    raise RuntimeError("Assertions must remain enabled; do not use -O or -OO.")

n, z, rho, u = s.symbols("n z rho u")
R = s.Rational
P = s.Matrix([
    [z*(4*n+1)**3*(16*n+13),
     -R(3,2)*z*(4*n+1)**3*(28*n**2+56*n+27),
     R(1,4)*z*(4*n+1)**3*(6*n+7)*(72*n**2+146*n+71)],
    [12*z*(96*n**3+132*n**2+50*n+5)-2*(576*n**3+876*n**2+404*n+59),
     -2*z*(1504*n**4+3856*n**3+3246*n**2+1018*n+97)
     +3024*n**4+8208*n**3+7632*n**2+2876*n+383,
     z*(7680*n**5+28864*n**4+40096*n**3+25020*n**2+6742*n+607)
     -R(1,2)*(6*n+7)*(2592*n**4+7128*n**3+6708*n**2+2570*n+353)],
    [4*(z-1)*(84*n**2+104*n+29),
     -2*(z-1)*(432*n**3+1044*n**2+776*n+173),
     z*(2176*n**4+7840*n**3+10128*n**2+5472*n+1017)
     -(6*n+7)*(360*n**3+876*n**2+650*n+143)]])
scale = 9*(2*n+1)**2*(6*n+5)**2*(6*n+7)**2/(16*z**3*(n+1)**4*(4*n+1)**3*(4*n+3)**3)
T = scale*P
expected_det = 729*(2*n+1)**6*(6*n+5)**6*(6*n+7)**6/(64*z**6*(n+1)**6*(4*n+1)**6*(4*n+3)**6)
assert s.cancel(scale**3*P.det()-expected_det) == 0
assert all(s.Poly(entry,z).degree() <= 1 for entry in P)

limit = s.Matrix(3,3,lambda i,j:s.limit((n+1)**i*T[i,j]/(n+2)**j,n,s.oo))
expected_limit = s.Matrix([
    [729/z**2, -15309/(8*z**2),19683/(4*z**2)],
    [6561*(z-1)/(8*z**3),-729*(188*z-189)/(64*z**3),2187*(80*z-81)/(32*z**3)],
    [15309*(z-1)/(64*z**3),-19683*(z-1)/(32*z**3),729*(136*z-135)/(64*z**3)]])
assert (limit-expected_limit).applyfunc(s.cancel) == s.zeros(3)
chi = rho**3-(3*z+R(27,2))*rho**2+(3*z**2-R(81,2)*z+R(729,16))*rho-z**3
assert s.cancel((limit-729/(16*z**3)*rho*s.eye(3)).det()+R(387420489,4096)/z**9*chi)==0
disc = s.factor(s.discriminant(chi,rho))
assert s.expand(disc-R(19683,256)*z*(1-z)*(32*z-81)**2)==0
assert s.expand(chi.subs(rho,2)-(R(1,2)-z)*(z**2-R(11,2)*z+R(361,4)))==0
e0=s.Matrix([1,0,0]); dual=limit.inv().T
cyclic=s.factor(s.Matrix.hstack(e0,dual*e0,dual**2*e0).det())
assert cyclic == -R(8192,14348907)*z**5*(32*z-81)

# Root bracketing lemma: all coefficients and comparisons are exact.
t1=R(16,729)*z**3*(1+R(8,9)*z+z**2)
residual=s.Poly(s.cancel((chi.subs(rho,t1)/z**5-R(67,243))/z),z)
assert residual.degree()==9 and all(abs(c)<1 for c in residual.all_coeffs())
zmax=R(1,7165); umax=R(119,10000)
assert zmax<umax**2 and R(67,243)-zmax/(1-zmax)>0
c=R('7.794229'); m2=R(27,4)-c*u+R(3,2)*u**2
coeffs=[R(154019907,4*10**12),-R(44461719470989,10**18),
        R(10125017113323,2*10**12),-R(23382687,4*10**6)]
target=u**2*(coeffs[0]+coeffs[1]*u+u**2*(coeffs[2]+coeffs[3]*u)+R(1,8)*u**4)
assert s.expand(chi.subs({z:u**2,rho:m2})-target)==0
assert coeffs[0]+coeffs[1]*umax>0 and coeffs[2]+coeffs[3]*umax>0
assert c>3*umax  # the middle-root lower bound increases with M
assert s.expand(chi.subs(rho,R(27,4))+z*(R(6561,16)-R(81,4)*z+z**2))==0

ctx.dps=80
M=arb(7165)
middle=arb(27)/4-(arb(7794229)/10**6)/M.sqrt()+arb(3)/(2*M)
margin=M.log()-arb(16).log()-8+middle.log()-(1+arb(8)/(9*M)+1/M**2).log()
assert margin>arb(89)/10**6
report={
    "status":"passed", "paper_labels":["eq:hh-T","eq:hh-detT","eq:hh-Tinf","eq:hh-cyclic","lem:hh-roots","prop:hh-gap"],
    "trajectory_determinant_identity":True,"limiting_matrix_identity":True,
    "characteristic_polynomial_identity":True,"discriminant":str(disc),
    "cyclic_determinant":str(cyclic),"root_bound_polynomial_identities":True,
    "small_root_residual_degree":residual.degree(),
    "small_root_residual_coefficients":[str(x) for x in residual.all_coeffs()],
    "middle_root_expansion_coefficients":[str(x) for x in coeffs],
    "threshold_M":7165,"positive_margin_arb":str(margin),
    "positive_margin_lower_rational":"89/1000000",
    "scope":"Exact algebra and a rigorous threshold enclosure, conditional on the manuscript's analytic and arithmetic arguments; no experimental figure data are generated."}
target=Path(__file__).with_name('gregory_certificate.json')
target.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
