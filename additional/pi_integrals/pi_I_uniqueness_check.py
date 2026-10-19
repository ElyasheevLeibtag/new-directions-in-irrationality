from pathlib import Path
import sympy as s,json
ns={};exec(Path(__file__).with_name('pi_I_axis_original.py').read_text(),ns)
A,B,G=ns['A'],ns['B'],ns['G'];y=s.symbols('y');H=y*y+6*y+25;v=y-25
# Coefficient map for T[H P/v^3] + m0 + m1*y + m2/v^2.
def T(P):return s.Poly(s.expand((1+2*G)*H*P*v+2*y*((s.diff(H,y)*P+H*s.diff(P,y))*v-3*H*P)+4*A*y*s.diff(H,y)*P*v-(4*B+20)*y*H*P),y)
polys=[s.Poly(v**4,y),s.Poly(y*v**4,y),s.Poly(v**2,y)]+[T(y**k) for k in range(7)]
mat=s.Matrix([[p.nth(k) for p in polys] for k in range(10)])
point={A:1,B:2,G:1};det=mat.subs(point).det();print('normal form determinant at seed',det)
out={'normal_form_determinant_at_I_seed':str(det),'normal_form_unique_generically':det!=0,'axis_determinants_at_I_seed':{name:str(ns[name].subs(point).det()) for name in ['MA','MB','MG']}}
print(out);Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))

assert det != 0 and all(ns[name].subs(point).det()!=0 for name in ["MA","MB","MG"]), "Seed is degenerate."
