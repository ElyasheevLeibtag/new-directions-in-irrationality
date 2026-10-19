from pathlib import Path
import sympy as s,json
m=s.symbols('m');Z={j:s.symbols('Z'+str(j)) for j in range(1,7)}
f=11664*(m+1)*(m+4)/((m+2)**6*(m+3)**6)
state=[]
for j in range(7):
 term=s.apart(m**j*f,m);total=0
 for summand in s.Add.make_args(term):
  num,den=s.fraction(summand);pol=s.Poly(den,m); roots=s.roots(den,m)
  assert len(roots)==1
  pole,k=next(iter(roots.items()));c=-pole;coeff=s.limit(summand*(m+c)**k,m,pole)
  total+=coeff*(Z[k]-s.harmonic(c-1,k))
 state.append(s.expand(total))
P=s.Matrix([[s.Rational(16342893984375,8),s.Rational(7888732953125,2),s.Rational(24573697734375,8),1206481062500,s.Rational(1912606640625,8),s.Rational(76504265625,4),0]])
Q=s.Matrix([[1969765000000,s.Rational(15220831421875,4),s.Rational(11858161171875,4),s.Rational(2329422984375,2),s.Rational(461663671875,2),18466546875,0]])
a=s.expand((P*s.Matrix(state))[0]);b=s.expand((Q*s.Matrix(state))[0]);print('Pstate',a,'Qstate',b,'difference',s.expand(a-Z[5]*b),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'seed':'5,2^7;4^7','state':[str(v) for v in state],'P':[str(v) for v in P],'Q':[str(v) for v in Q],'Pstate':str(a),'Qstate':str(b),'identity':s.expand(a-Z[5]*b)==0},indent=2))

assert s.expand(a-Z[5]*b)==0
