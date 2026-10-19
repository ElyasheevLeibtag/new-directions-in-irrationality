import sympy as s, json, time
from pathlib import Path
ns={};exec(Path(__file__).with_name('pi_I_axis_original.py').read_text(),ns)
A,B,G=ns['A'],ns['B'],ns['G']; y=s.symbols('y');H=y*y+6*y+25;v=y-25
basis=[1,y,1/v**2];factors=[H**2,1/v**2,y]
Mlist=[ns['MA'],ns['MB'],ns['MG']]
def Tnum(P):
 return s.Poly(s.expand((1+2*G)*H*P*v+2*y*((s.diff(H,y)*P+H*s.diff(P,y))*v-3*H*P)+4*A*y*s.diff(H,y)*P*v-(4*B+20)*y*H*P),y)
ops=[Tnum(y**k) for k in range(7)]
checks=[];certs=[]
for ai,(F,M) in enumerate(zip(factors,Mlist)):
 for j in range(3):
  tt=time.time();target=s.Poly(s.cancel((F*basis[j]-sum(M[k,j]*basis[k] for k in range(3)))*v**4),y)
  # triangular coefficient reduction in Q(A,B,G)[y]
  coeff={i:s.cancel(target.nth(i)) for i in range(10)};p={}
  for k in reversed(range(7)):
   val=s.cancel(coeff[k+3]/ops[k].nth(k+3));p[k]=val
   for i in range(k+4): coeff[i]=s.cancel(coeff[i]-val*ops[k].nth(i))
  rem=[s.cancel(coeff[i]) for i in range(10)]
  ok=all(z==0 for z in rem);print(ai,j,ok,round(time.time()-tt,2),flush=True)
  checks.append(ok);certs.append({'axis':str([A,B,G][ai]),'column':j,'P_coefficients':[str(p[k]) for k in range(7)],'verified':ok,'remainder':[str(z) for z in rem]})
Path(__file__).with_suffix('.json').write_text(json.dumps({'operator':'T[h]=h+2y h_prime+2y h*(2A H_prime/H+G/y-(2B+10)/(y-25))','H':'y^2+6y+25','h':'H*P/(y-25)^3','all_nine_identities':all(checks),'certificates':certs},indent=2))

assert all(checks), "A symbolic axis identity failed."
