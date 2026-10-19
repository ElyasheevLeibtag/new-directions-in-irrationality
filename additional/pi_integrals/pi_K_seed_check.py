import sympy as s,json
from pathlib import Path
x=s.symbols('x');H=((x-4)**2+4)*((x-6)**2+4);lo=4-2*s.I;hi=4+2*s.I

def pair(a,b,g):
 num=s.Poly(H**a*(x-5)**(2*g),x);den=s.Poly(x**b*(x-10)**b,x)
 quo,rem=s.div(num,den);f=rem.as_expr()/den.as_expr();c={}
 for v in (0,10):
  y=s.symbols('y'); series=s.Poly(s.series(s.cancel(f*(x-v)**b).subs(x,v+y),y,0,b).removeO(),y)
  for k in range(1,b+1):c[v,k]=series.nth(b-k)
 assert c[0,1]==-c[10,1]
 C=c[0,1]/2
 R=s.integrate(quo.as_expr(),x)+sum(c[v,k]*(x-v)**(1-k)/(1-k) for v in (0,10) for k in range(2,b+1))
 D=s.simplify(-s.I*(R.subs(x,hi)-R.subs(x,lo)))
 return C,D
pairs=[pair(0,1,0),pair(0,1,1),pair(0,3,0)]
print(pairs,flush=True)
P=s.Matrix([[-z[1] for z in pairs]]);Q=s.Matrix([[z[0] for z in pairs]])
L=s.ilcm(*[s.denom(v) for v in P],*[s.denom(v) for v in Q]);P=L*P;Q=L*Q;gc=s.igcd(*[int(z) for z in P],*[int(z) for z in Q]);P=P/gc;Q=Q/gc
print('P',P,'Q',Q,flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'a_b_g':[0,1,0],'basis':'K(a,b,g),K(a,b,g+1),K(a,b+2,g)','pairs':[[str(a),str(b)] for a,b in pairs],'P':[str(x) for x in P],'Q':[str(x) for x in Q]},indent=2))
