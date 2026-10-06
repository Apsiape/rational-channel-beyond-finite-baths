#!/usr/bin/env python3
"""Exact checks for Appendix B (root SOS identity, Netzer--Thom transfer constant) and for the
final exponent comparison in Lemma C.1 (8 M0 < J), plus the 12-qubit counts.

Standard library only.  It checks scalar identities and inequalities; it does not prove the
analytic statements, which are Wang and Zhi's or are proved in the paper.
"""
from fractions import Fraction as F


def check(cond,msg):
    if not cond: raise RuntimeError(msg)

# ---------- tiny polynomial algebra over Q ----------
def trim(p):
    p=list(p)
    while p and p[-1]==0:p.pop()
    return p or [F(0)]
def add(a,b):
    n=max(len(a),len(b));r=[F(0)]*n
    for i in range(n):r[i]=(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
    return trim(r)
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    r=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return trim(r)
def divmodp(a,b):
    a=trim(a);b=trim(b);q=[F(0)]*max(1,len(a)-len(b)+1)
    while not (len(a)==1 and a[0]==0) and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]/b[-1];q[k]+=c
        sub=[F(0)]*k+scale(b,c);a=add(a,scale(sub,-1))
    return trim(q),trim(a)
T=[F(0),F(1)]; one=[F(1)]
def sub(a,b):return add(a,scale(b,-1))
def powp(a,n):
    r=one
    while n:
        if n&1:r=mul(r,a)
        a=mul(a,a);n//=2
    return r

t1=sub(T,one)
E1=mul(mul(T,sub(T,[F(2)])), sub(scale(powp(t1,2),8),one))
E2=scale(mul(mul(T,t1), sub(scale(powp(t1,2),8),one)),F(1,14))
Estar=scale(mul(mul(T,sub(T,[F(2)])),powp(t1,2)),F(-64,7))
Vstar=mul(sub(T,[F(21,32)]),Estar)
lhs=sub(powp(T,2),scale(T,F(5,8)))
rhs=add(add(scale(powp(E1,2),F(3,8)),scale(powp(E2,2),F(11,4))),add(scale(powp(Vstar,2),F(2)),scale(powp(Estar,2),F(7,512))))
m=mul(mul(mul(T,t1),sub(T,[F(2)])),sub(powp(t1,2),[F(1,8)]))
q,rem=divmodp(sub(lhs,rhs),m)
check(rem==[F(0)],'finite-field spectral polynomial identity failed')

# ---------- reconstructed integer certificate headline constants ----------
nt_l1=F(572700482439,25000000000000)
nt_eff=F(561,2000)-4*nt_l1
check(nt_eff==F(1180424517561,6250000000000),'NT effective gap mismatch')
check(nt_eff>F(1,6),'NT corrected gap is not >1/6')
CK=F(37610695003407823,28800000000000000)
transfer=14*CK
check(transfer==F(263274865023854761,14400000000000000),'NT transfer mismatch')
check(transfer<19,'NT transfer not <19')

# ---------- WZ coarse epsilon_in comparison, checked at exponent level ----------
# T0 = 2^1121200; M0 = 2^(32768*T0); epsilon_in >= 2^(-8M0).
# J = 2^(2^1122000), delta=2^-J.  It suffices that 8M0 < J,
# i.e. 3 + 32768*T0 < 2^1122000.
T0_exp=1121200
lhs_exp=3 + (1 << (T0_exp+15))
rhs_exp=(1 << 1122000)
check(lhs_exp<rhs_exp,'WZ exponent comparison 8M0<J failed')
# Smaller coarse inequalities stated in the reconstruction.
check(1121072 + 8*(1<<1121111) < 16*(1<<1121200),'gamma last-entry exponent bound failed')
check(16384*(1<<1121200)+10 <= 2*(1 << (15+1121200)),'mu_I/M0 exponent bound failed')

# ---------- dimension/count arithmetic and exponent formulas ----------
labels=949; triangles=1866; d12=4096; anchors=d12-2*triangles
check(anchors==364,'12q anchor count')
check(4*triangles+anchors==7828,'12q nonzero count')

print('OK exact rigidity/compiler scalar checks')
print('finite-field SOS polynomial identity: exact')
print('NT effective gap =',nt_eff,'>',F(1,6))
print('NT relation transfer =',transfer,'< 19')
print('12q anchors=364 labels/rank=949 triangles=1866 nonzero=7828')
print('WZ tolerance exponent comparison: 8 M0 < J (Lemma C.1)')
