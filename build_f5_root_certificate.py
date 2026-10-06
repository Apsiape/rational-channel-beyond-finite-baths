from fractions import Fraction
import sympy as sp, math, json, os
x=sp.symbols('x')
Q=sp.QQ
c=sp.Rational(11,20)
# irreducible/rational factors for spectra r=0..4
factors=[x, x-1, x-2, x-sp.Rational(6,5), x-sp.Rational(4,5), x-sp.Rational(26,25), x-sp.Rational(24,25),
         (x-1)**2-sp.Rational(1,5), (x-1)**2-sp.Rational(1,125)]
m=sp.Poly(1,x,domain=Q)
for f in factors: m*=sp.Poly(f,x,domain=Q)
m=sp.monic(m)
q=sp.Poly(x*(x-c),x,domain=Q)

def inv_mod_poly(a,f):
    s,t,g=sp.gcdex(sp.Poly(a,x,domain=Q), sp.Poly(f,x,domain=Q))
    assert g.as_expr()==1
    return sp.rem(s, sp.Poly(f,x,domain=Q), domain=Q)

def crt_idempotent(f):
    f=sp.Poly(f,x,domain=Q)
    Mi=sp.div(m,f)[0]
    inv=inv_mod_poly(Mi,f)
    e=sp.rem(Mi*inv,m,domain=Q)
    # verify
    assert sp.rem(e-1,f,domain=Q).is_zero
    assert sp.rem(e*e-e,m,domain=Q).is_zero
    return e

def four_squares(n):
    n=int(n); lim=math.isqrt(n)
    pair={}
    for i in range(lim+1):
        ii=i*i
        for j in range(i,lim+1):
            s=ii+j*j
            if s<=n and s not in pair: pair[s]=(i,j)
    for s,ij in pair.items():
        if n-s in pair:
            return (*ij,*pair[n-s])
    raise RuntimeError(n)

def rational_squares(fr):
    fr=Fraction(fr)
    N=fr.numerator*fr.denominator
    vals=four_squares(N)
    D=fr.denominator
    return [sp.Rational(v,D) for v in vals if v]

squares=[]
entries=[]
for idx,fexpr in enumerate(factors):
    f=sp.Poly(fexpr,x,domain=Q)
    e=crt_idempotent(f)
    deg=f.degree()
    if f.as_expr()==x:
        reps=[] # q=0 at x=0
    elif deg==1:
        root=sp.solve(f.as_expr(),x)[0]
        val=sp.Rational(sp.simplify(q.as_expr().subs(x,root)))
        reps=[sp.Poly(r,x,domain=Q) for r in rational_squares(Fraction(int(val.p),int(val.q)))]
    else:
        # explicit Q(sqrt5) decompositions described in proof
        if sp.expand(f.as_expr())==sp.expand((x-1)**2-sp.Rational(1,5)):
            # q = (29/50 + sqrt5/4)^2 + (3/100)^2+(1/100)^2+(1/100)^2; sqrt5=5(x-1)
            reps=[sp.Poly(sp.Rational(29,50)+sp.Rational(5,4)*(x-1),x,domain=Q),
                  sp.Poly(sp.Rational(3,100),x,domain=Q),sp.Poly(sp.Rational(1,100),x,domain=Q),sp.Poly(sp.Rational(1,100),x,domain=Q)]
        else:
            # q = (29/100 + sqrt5/10)^2 + 1^2+3^2+27^2+50^2 over 100^2; sqrt5=25(x-1)
            reps=[sp.Poly(sp.Rational(29,100)+sp.Rational(5,2)*(x-1),x,domain=Q),
                  sp.Poly(sp.Rational(1,100),x,domain=Q),sp.Poly(sp.Rational(3,100),x,domain=Q),
                  sp.Poly(sp.Rational(27,100),x,domain=Q),sp.Poly(sp.Rational(50,100),x,domain=Q)]
        chk=sp.rem(q - sum((a*a for a in reps),sp.Poly(0,x,domain=Q)), f, domain=Q)
        assert chk.is_zero, chk
    lifted=[]
    for a in reps:
        p=sp.rem(e*a,m,domain=Q)
        squares.append(p); lifted.append(p)
    entries.append({'factor':str(f.as_expr()),'idempotent':str(e.as_expr()),'n_squares':len(lifted)})

S=sum((p*p for p in squares),sp.Poly(0,x,domain=Q))
rem=sp.rem(q-S,m,domain=Q)
assert rem.is_zero, rem
h=sp.div(q-S,m)[0]

def l1(poly):
    return sum(abs(sp.Rational(c)) for c in poly.all_coeffs())
print('m degree',m.degree(),'squares',len(squares),'identity exact')
print('m=',sp.factor(m.as_expr()))
print('h degree',h.degree(),'h l1',l1(h))
print('m l1',l1(m),'q l1',l1(q))
print('max square l1',max(l1(p) for p in squares),'sum square l1sq',sum(l1(p)**2 for p in squares))
# multiplicativity residual constant for m(A) using ||A||_1<=4
coeffs={i:sp.Rational(m.nth(i)) for i in range(m.degree()+1)}
Cm=sum(abs(v)*max(0,i-1)*4**i for i,v in coeffs.items())
print('Cm',Cm, float(Cm))
# functional h(A) op norm bound from coeff l1 * ||A||^k, ||A||<=4
Ch=sum(abs(sp.Rational(h.nth(i)))*4**i for i in range(h.degree()+1))
print('Ch_op_bound',Ch,float(Ch))
# save compact
out={'c':[11,20],'m_coeff_low_to_high':[[int(sp.numer(m.nth(i))),int(sp.denom(m.nth(i)))] for i in range(m.degree()+1)],
     'h_coeff_low_to_high':[[int(sp.numer(h.nth(i))),int(sp.denom(h.nth(i)))] for i in range(h.degree()+1)],
     'squares':[],'Cm':[int(sp.numer(Cm)),int(sp.denom(Cm))],'Ch':[int(sp.numer(Ch)),int(sp.denom(Ch))], 'entries':entries}
for p in squares:
    out['squares'].append([[int(sp.numer(p.nth(i))),int(sp.denom(p.nth(i)))] for i in range(m.degree())])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'f5_root_sos.json'),'w') as f: json.dump(out,f,separators=(',',':'),sort_keys=True)
