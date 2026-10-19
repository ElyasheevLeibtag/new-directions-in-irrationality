import sympy as s, json
from pathlib import Path
x=s.symbols('x');ii=s.I
z0=4-2*ii;z1=4+2*ii
H=((x-4)**2+4)*((x-6)**2+4)
def pair(A,B,G):
 f=H**(2*A)*(x-5)**(2*G)/(x**(2*B+10)*(x-10)**(2*B+10))
 f=s.cancel(f);b=2*B+10
 # Laurent coefficient at each pole from truncated Taylor expansion.
 c={}
 for v in (0,10):
  y=s.symbols('y'); series=s.Poly(s.series(s.cancel(f*(x-v)**b).subs(x,v+y),y,0,b).removeO(),y)
  for k in range(1,b+1): c[v,k]=series.nth(b-k)
 assert c[0,1]==-c[10,1]
 C=c[0,1]/2
 # Rational antiderivative; there is no polynomial part for these initial values.
 R=sum(c[v,k]*(x-v)**(1-k)/(1-k) for v in (0,10) for k in range(2,b+1))
 D=s.simplify(-ii*(R.subs(x,z1)-R.subs(x,z0)))
 return C,D
pairs=[pair(1,2,1),pair(1,2,2),pair(1,3,1)]
for p in pairs: print(p,flush=True)
P=s.Matrix([[-5236454140000,249191746500000,3609168431]])
Q=s.Matrix([[-4799995200000,83002920000000,-3414290880]])
C=s.Matrix([[v[0] for v in pairs]]);D=s.Matrix([[v[1] for v in pairs]])
scale=s.factor(C[0]/Q[0]); print('scale',scale);print('C-Qscale',s.simplify(C-scale*Q));print('D+Pscale',s.simplify(D+scale*P))
Path(__file__).with_suffix('.json').write_text(json.dumps({'basis_pairs':[[str(a),str(b)] for a,b in pairs],'scale':str(scale),'coefficient_identity':C==scale*Q and D==-scale*P},indent=2))

assert C==scale*Q and D==-scale*P, "Initial-state identity failed."
