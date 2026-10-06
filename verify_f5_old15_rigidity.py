#!/usr/bin/env python3
"""Exact scalar/polynomial checks for the characteristic-five quantitative rigidity route.
Standard library only.  Analytic lemmas (pair spectrum and collection algorithm) are proved in the paper;
this verifier checks the shipped SOS identity, its coefficient ledger, the 13-cell Steinberg reductions,
and all numerical slack used to turn them into an explicit old-15 tolerance.
"""
from fractions import Fraction as F
import json, os
HERE=os.path.dirname(__file__)

def ck(c,m):
    if not c: raise RuntimeError(m)

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

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
def sub(a,b):return add(a,scale(b,-1))
def powp(a,n):
    r=[F(1)]
    while n:
        if n&1:r=mul(r,a)
        a=mul(a,a);n//=2
    return r

def pair(a,b):return F(int(a),int(b))
D=json.load(open(os.path.join(HERE,'f5_root_sos.json')))
m=[pair(*z) for z in D['m_coeff_low_to_high']]
h=[pair(*z) for z in D['h_coeff_low_to_high']]
squares=[[pair(*z) for z in p] for p in D['squares']]
# Intended annihilating polynomial t(t-1)(t-2)(t-6/5)(t-4/5)(t-26/25)(t-24/25)
# *((t-1)^2-1/5)*((t-1)^2-1/125)
t=[F(0),F(1)];one=[F(1)]
mi=[F(1)]
for root in [F(0),F(1),F(2),F(6,5),F(4,5),F(26,25),F(24,25)]:
    mi=mul(mi,[-root,F(1)])
mi=mul(mi,sub(powp(sub(t,one),2),[F(1,5)]))
mi=mul(mi,sub(powp(sub(t,one),2),[F(1,125)]))
ck(trim(mi)==trim(m),'F5 annihilating polynomial mismatch')
q=sub(powp(t,2),scale(t,F(11,20)))
S=[F(0)]
for p in squares:S=add(S,mul(p,p))
ck(sub(q,S)==mul(h,m),'F5 SOS polynomial identity failed')
ck(len(squares)==23,'unexpected F5 square count')
# coefficient ledgers used in approximate transfer
Cm=sum(abs(c)*max(0,i-1)*(4**i) for i,c in enumerate(m))
Lm=sum(abs(c)*i*(4**(i-1)) for i,c in enumerate(m) if i)
Ch=sum(abs(c)*(4**i) for i,c in enumerate(h))
ck([Cm.numerator,Cm.denominator]==D['Cm'],'Cm mismatch')
ck([Ch.numerator,Ch.denominator]==D['Ch'],'Ch mismatch')
A=1<<20 # proved collection-area ceiling per product of canonical pair words
Croot=F(2*A,36)*(3*Ch*(Cm+Lm)+9)
ck(Croot < (1<<114),'root residual constant is not <2^114')
# ---- verify 3 opposite-root relations from the first 12 Steinberg balanced equations ----
def fr(w):
    st=[]
    for z in w:
        if st and st[-1]==-z:st.pop()
        else:st.append(z)
    return tuple(st)
def invw(w):return tuple(-z for z in reversed(w))
def rotations(w):
    for k in range(len(w)):yield w[k:]+w[:k]
def cyc_reduce(w):
    w=fr(w)
    while len(w)>=2 and w[0]==-w[-1]:w=fr(w[1:-1])
    return w
def cyc_can(w):
    w=cyc_reduce(w)
    if not w:return ()
    vals=[]
    for v in (w,invw(w)):vals.extend(rotations(v))
    return min(vals,key=lambda z:(len(z),z))
rels=[
 ((1,2),(2,1)),((3,4),(4,3)),((5,6),(6,5)),
 ((3,5),(5,3)),((1,6),(6,1)),((2,4),(4,2)),
 ((1,4),(2,4,1)),((2,6),(1,6,2)),((3,2),(4,2,3)),
 ((4,5),(3,5,4)),((5,1),(6,1,5)),((6,3),(5,3,6)),
]
# ---- the twelve relations above are the earlier report's copy-zero Steinberg relators ----
# Script index k (1..6) is U_ab in the order (1,2),(1,3),(2,1),(2,3),(3,1),(3,2), which is
# the earlier report's t_12,t_13,t_21,t_23,t_31,t_32 of copy 0, generator ids 22..27.
P=json.load(open(os.path.join(HERE,'earlier','presentation.json')))
names={g['id']:g['name'] for g in P['generators']}
def to_v1(w):return tuple((21+abs(z))*(1 if z>0 else -1) for z in w)
ours={cyc_can(to_v1(fr(tuple(l)+invw(tuple(r))))) for l,r in rels}
steinberg=[tuple(r['word']) for r in P['relators'] if r['family'].startswith('action-Steinberg')]
copy0={cyc_can(w) for w in steinberg if all(22<=abs(z)<=27 for z in w)}
ck(len(ours)==12 and ours==copy0,'Steinberg relations differ from the earlier presentation, copy 0')
from collections import Counter
fam=Counter(r['family'] for r in P['relators'])
ck(dict(fam)=={'triangle-order':30,'triangle-abelian':108,'triangle-class-two':3732,
               'action-Steinberg-commute':12,'action-Steinberg-product':12,
               'action-conjugation':504,'negative-name':18,'witness-name':1},f'relator families {dict(fam)}')
ck(P['witness_word']==[3,37,-22,-3,22,-37] and P['witness_generator']==43,'witness word')
ck(any(r['family']=='witness-name' and r['word']==[-43,3,37,-22,-3,22,-37] for r in P['relators']),'witness relator')
ck(names[3]=='a2','h = a_2')
# forward compression relators t g t^-1 = w_g for the six action generators of each copy and the
# twelve positive roots g: present, targets are positive-root words of length 1 or 16, and the two
# copies (ids 22..27 and 37..42, in the same order) have identical targets (Theorem D.4).
conj=[tuple(r['word']) for r in P['relators'] if r['family']=='action-conjugation']
targets={}
for w in conj:
    t_,g_=w[0],w[1]
    if t_>0 and len(w)>=4 and w[2]==-t_ and 1<=g_<=12 and (22<=t_<=27 or 37<=t_<=42):
        tgt=invw(w[3:])
        ck(all(1<=abs(z)<=12 for z in tgt) and len(tgt) in (1,16),f'compression target {w}')
        ck((t_,g_) not in targets,'duplicate compression relator')
        targets[(t_,g_)]=tgt
ck(len(targets)==144,f'forward compression relators {len(targets)}')
ck(all(targets[(t_,g_)]==targets[(t_+15,g_)] for t_ in range(22,28) for g_ in range(1,13)),'copies differ')
rules=set()
for rid,(lhs,rhs) in enumerate(rels):
    R=fr(tuple(lhs)+invw(tuple(rhs)))
    for ori in (R,invw(R)):
        for rot in rotations(ori):
            for k in range(1,len(rot)+1):
                a=rot[:k]; b=invw(rot[k:]); rules.add((rid,a,b))
def can_transition(before,after,rid,a,b):
    before=tuple(before);after=tuple(after);a=tuple(a);b=tuple(b)
    if (rid,a,b) not in rules:return False
    for v0 in (before,invw(before)):
        for v in rotations(v0):
            for i in range(len(v)-len(a)+1):
                if v[i:i+len(a)]==a:
                    if cyc_can(fr(v[:i]+b+v[i+len(a):]))==after:return True
    return False
R=json.load(open(os.path.join(HERE,'st3_opposite_reductions.json')))
ck(len(R['records'])==3,'need 3 opposite reductions')
for rec in R['records']:
    cur=cyc_can(tuple(rec['target']))
    ck(tuple(rec['steps'][0][0])==cur,'bad reduction start')
    for step in rec['steps']:
        before,after,rid,flag,a,b=step
        ck(tuple(before)==cur,'noncontiguous reduction')
        ck(can_transition(cur,tuple(after),rid,a,b),'invalid Steinberg cell')
        cur=tuple(after)
    ck(cur==() and rec['cost']==len(rec['steps'])==13,'bad reduction end/cost')
# integer residual coefficient: reconstructed NT transfer <19 times derived-relator cost 13
Cint=19*13
ck(Cint==247,'integer coefficient')
# Same colossal delta=2^{-J} lies below epsilon_in for root and integer spectral inputs.
# WZ exponent proof gives eps_in >= 2^{-8 M0} for h<=2^14, 2^-9<=kappa<=1/4.
# C only enters as /(C+1); C<=2^114 still fits because 114 << 4 M0.
T0=1<<1121200
M0=1<<(15+1121200)  # exponent of 2 in M0: M0 itself is 2^(32768*T0); we compare exponents symbolically below
# avoid constructing 2^(32768*T0): same comparison as base verifier
ck(3 + (1 << (1121200+15)) < (1 << 1122000),'8M0<J exponent comparison')
ck(115 < 4*(1<<15),'C absorption sanity') # vastly weaker than true 115<4M0
# Fixed reconstruction constants
k=F(1,1<<9); r_exp=17335; eta_exp=17313
# positive-root averaging: 3*5^4=1875<2^11; lambda=1/60 gives kappa=1/120 between 2^-9 and 1/4
ck(3*5**4 < (1<<11),'root averaging length')
ck(F(1,1<<9) <= F(1,120) <= F(1,4),'root kappa range')
ck(Croot < (1<<114),'root C bound')
ck(Cint < (1<<8),'integer C bound')
# Near inclusion ledger: base-root commutator <=r; compression target length <=16; e<=delta<<r.
# generator commutator <=18r, subgroup word length<=16 => energy root <=288r.
# k^{-1/2}<32, hence distance <(32*288+1)r <2^14 r <2^22 r=eta.
ck(32*288+1 < (1<<14),'near inclusion coefficient')
# eta/r = 2^(17335-17313)=2^22
ck((1<<14) < (1<<22),'near inclusion fits reversal eta')
# common compressor: <=2e on 12 base roots; subgroup word <=16 => energy root <=32e;
# Poincare <=32*32 e + r <2^11 r (as e<<r).
ck(32*32+1 < (1<<11),'compressor-to-A coefficient')
# reverse inclusion adds 2^-20; commutator with h <= r + 2*(2^11 r +2^-20) < 2^-18.
# verify using rational powers at worst r=2^-17335.
r=F(1,1<<r_exp)
witness = r + 2*((1<<11)*r + F(1,1<<20))
ck(witness < F(1,1<<18),'witness bound')
# named witness relation adds e<delta<2^-18, still <1/2.
print('OK characteristic-five rigidity scalars and certificates')
print('F5 root SOS: degree',len(m)-1,'squares',len(squares))
print('C_root =',Croot,'< 2^114')
print('Steinberg opposite reductions: 3 x 13 cells exact')
print('integer residual coefficient <=247')
print('Steinberg relations = earlier presentation copy-zero relators; relator families and witness word match')
print('144 forward compression relators: targets of length 1 or 16 in positive roots, identical in both copies')
print('43-generator presentation tolerance delta_5 = 2^-J, J = 2^(2^1122000) (Theorem D.4)')
