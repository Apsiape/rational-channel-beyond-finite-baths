#!/usr/bin/env python3
"""Independent exact verifier for the reconstructed Netzer--Thom spectral certificate.
Uses only Python's standard library.  It verifies:
  * the rational SOS Gram factor Q and augmentation normalization;
  * the exact group-algebra residual in SL(3,Z);
  * the Netzer--Thom correction margin > 1/6;
  * the exact residual coefficient-transfer constant;
  * all finite word-reduction certificates from the 15 balanced relations.
"""
from fractions import Fraction
from collections import defaultdict, Counter
import json, hashlib, os, sys
from datafile import read_data

ROOT=os.path.dirname(os.path.abspath(__file__))
SOS=os.path.join(ROOT,'nt_sos_certificate.json.gz')
RED=os.path.join(ROOT,'nt_word_reductions.json.gz')

I3=((1,0,0),(0,1,0),(0,0,1))
def mmul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
def minv(A):
    cof=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rs=[r for r in range(3) if r!=i]; cs=[c for c in range(3) if c!=j]
            d=A[rs[0]][cs[0]]*A[rs[1]][cs[1]]-A[rs[0]][cs[1]]*A[rs[1]][cs[0]]
            if (i+j)&1:d=-d
            row.append(d)
        cof.append(row)
    return tuple(tuple(cof[j][i] for j in range(3)) for i in range(3))

def fr(w):
    st=[]
    for x in w:
        if st and st[-1]==-x: st.pop()
        else: st.append(x)
    return tuple(st)
def invw(w): return tuple(-x for x in reversed(w))
def rotations(w):
    for k in range(len(w)): yield w[k:]+w[:k]
def cyc_reduce(w):
    w=fr(w)
    while len(w)>=2 and w[0]==-w[-1]: w=fr(w[1:-1])
    return w
def cyc_can(w):
    w=cyc_reduce(w)
    if not w:return ()
    cands=list(rotations(w))+list(rotations(invw(w)))
    return min(cands,key=lambda z:(len(z),z))

def parse_key(s):
    if not s:return ()
    return tuple(int(x) for x in s.split(','))

# Six positive elementary transvections in WZ order 12,13,21,23,31,32.
GEN=[]
for i,j in [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]:
    A=[list(r) for r in I3]; A[i][j]+=1; GEN.append(tuple(tuple(r) for r in A))
def evalw(w):
    A=I3
    for s in w:
        g=GEN[s-1] if s>0 else minv(GEN[-s-1])
        A=mmul(A,g)
    return A

# Natural basis: identity, 12 signed generators, then new signed-generator products.
signed=[]
for idx,g in enumerate(GEN,1): signed.extend([(g,(idx,)),(minv(g),(-idx,))])
B=[I3]; W=[()]; seen={I3:0}
for g,w in signed:
    if g not in seen: seen[g]=len(B); B.append(g); W.append(w)
for g,w in signed:
    for h,v in signed:
        x=mmul(g,h); ww=fr(w+v)
        if x not in seen: seen[x]=len(B); B.append(x); W.append(ww)
assert len(B)==121
BI=[minv(x) for x in B]

sos=json.loads(read_data(SOS))
Q=sos['Q']; den=int(sos['denominator']); eps=Fraction(*sos['epsilon'])
assert len(Q)<=121 and all(len(r)==121 for r in Q)
assert all(sum(map(int,r))==0 for r in Q), 'every SOS row must have augmentation zero'

# Exact P = Q^T Q.
P=[[0]*121 for _ in range(121)]
for row in Q:
    nz=[(i,int(x)) for i,x in enumerate(row) if int(x)]
    for i,x in nz:
        Pi=P[i]
        for j,y in nz: Pi[j]+=x*y
DD=den*den

# b = a^* P a in the actual matrix group.
b_group=defaultdict(Fraction)
for i in range(121):
    for j in range(121):
        c=P[i][j]
        if c: b_group[mmul(BI[i],B[j])]+=Fraction(c,DD)

# a = Delta^2 - eps Delta in the group algebra.
Delta=defaultdict(Fraction); Delta[I3]=Fraction(12)
for g,w in signed: Delta[g]-=1
a_group=defaultdict(Fraction)
for g,c in Delta.items():
    for h,d in Delta.items(): a_group[mmul(g,h)]+=c*d
for g,c in Delta.items(): a_group[g]-=eps*c
c_group={g:a_group.get(g,Fraction())-b_group.get(g,Fraction()) for g in set(a_group)|set(b_group)}
c_group={g:c for g,c in c_group.items() if c}
l1=sum((abs(c) for c in c_group.values()),Fraction())
assert l1==Fraction(*sos['l1']), (l1,sos['l1'])
effective=eps-4*l1
assert effective>Fraction(1,6)

# Build the free residual R = A_F - B_F - C_F using one representative per matrix class.
b_free=defaultdict(Fraction)
for i in range(121):
    wi_inv=invw(W[i])
    for j in range(121):
        c=P[i][j]
        if c:b_free[fr(wi_inv+W[j])]+=Fraction(c,DD)
Dfree=defaultdict(Fraction);Dfree[()]=Fraction(12)
for g,w in signed:Dfree[w]-=1
Afree=defaultdict(Fraction)
for u,c in Dfree.items():
    for v,d in Dfree.items():Afree[fr(u+v)]+=c*d
for u,c in Dfree.items():Afree[u]-=eps*c
E=defaultdict(Fraction);E.update(Afree)
for u,c in b_free.items():E[u]-=c
E={u:c for u,c in E.items() if c}
classes=defaultdict(dict)
for u,c in E.items():classes[evalw(u)][u]=c
canonical=set();sum_nonrep=Fraction();res_l1=Fraction();comparison_pairs=0
for g,d in classes.items():
    rep=min(d,key=lambda w:(len(w),w)); total=sum(d.values(),Fraction())
    rd=dict(d);rd[rep]=rd.get(rep,Fraction())-total
    rd={u:c for u,c in rd.items() if c}
    res_l1+=sum((abs(c) for c in rd.values()),Fraction())
    for u,c in rd.items():
        if u!=rep:
            comparison_pairs+=1; sum_nonrep+=abs(c)
            cw=cyc_can(invw(rep)+u)
            if cw:canonical.add(cw)
assert comparison_pairs==6722
Cscaled=Fraction(2)*sum_nonrep/Fraction(24*24)
assert Cscaled < Fraction(131,100)  # < 1.31
assert Cscaled == Fraction(37610695003407823, 28800000000000000)  # C_K of Appendix B

# Load and verify word certificates.
red=json.loads(read_data(RED))
base_rel=[(tuple(x[0]),tuple(x[1])) for x in red['base_relations']]
assert len(base_rel)==15
# Full one-cell split rules, exactly as used by the certificates.
full_rules=[];rseen=set()
for rid,(lhs,rhs) in enumerate(base_rel):
    R=fr(lhs+invw(rhs))
    for ori in (R,invw(R)):
        for rot in rotations(ori):
            L=len(rot)
            for k in range(1,L+1):
                a=rot[:k]; b=invw(rot[k:]); key=(a,b)
                if key not in rseen:rseen.add(key);full_rules.append((rid,a,b))

def apply_base_action(x,act):
    rid,flip,k,i,aa,bb=act; rid=int(rid); flip=int(flip); k=int(k); i=int(i)
    aa=tuple(aa);bb=tuple(bb)
    v0=invw(x) if flip else x
    assert v0
    v=v0[k:]+v0[:k]
    assert v[i:i+len(aa)]==aa
    assert any(rr==rid and A==aa and BB==bb for rr,A,BB in full_rules)
    return cyc_can(fr(v[:i]+bb+v[i+len(aa):]))

base_paths={parse_key(k):v for k,v in red['base_identity_paths'].items()}
base_cost={}
for w,path in base_paths.items():
    x=cyc_can(w)
    assert x==w
    for act in path:x=apply_base_action(x,act)
    assert x==(),('bad base path',w,x)
    base_cost[w]=len(path)
assert max(base_cost.values())<=10

# Canonical base relator words for direct macro steps.
base_can=[cyc_can(fr(lhs+invw(rhs))) for lhs,rhs in base_rel]
comp_paths={parse_key(k):v for k,v in red['comparison_paths'].items()}
assert set(comp_paths)==canonical, (len(comp_paths),len(canonical),len(set(comp_paths)^canonical))
assert len(canonical)==858
maxcost=0
for start,rec in comp_paths.items():
    x=start; total=0
    for step in rec['steps']:
        before=tuple(step[0]); after=tuple(step[1]); cost=int(step[2]); action=step[3]
        assert before==x
        rot_i=int(action[0]); aa=tuple(action[1]); bb=tuple(action[2]); acost=int(action[3]); cert=action[4]; mid=tuple(action[5])
        assert cost==acost
        v=x[rot_i:]+x[:rot_i]
        assert v[:len(aa)]==aa
        got=cyc_can(fr(bb+v[len(aa):]))
        assert got==after,(start,x,got,after)
        assert cyc_can(fr(aa+invw(bb)))==mid
        if cert[0]=='base':
            rid=int(cert[1]); assert cost==1; assert mid==base_can[rid]
        elif cert[0]=='solved':
            mw=tuple(cert[1]); assert mw==mid; assert mw in base_cost; assert cost==base_cost[mw]
        else: raise AssertionError(cert)
        total+=cost; x=after
    assert x==()
    assert total==int(rec['cost'])
    maxcost=max(maxcost,total)
assert maxcost==14

# Spectral transfer constant: conjugation costs a factor 2.
transfer = Fraction(maxcost)*Cscaled
assert transfer < Fraction(19)  # in fact <18.3

print('OK exact Netzer--Thom reconstruction')
print('basis_size=121 sos_rows=%d denominator=%d' % (len(Q),den))
print('correction_l1=%s ~= %.12g' % (l1,float(l1)))
print('effective_Delta_gap=%s ~= %.12g (>1/6)' % (effective,float(effective)))
print('residual_comparisons=%d canonical=%d max_relation_replacements=%d' % (comparison_pairs,len(canonical),maxcost))
print('C_K=%s ~= %.12g; spectral_transfer=%s ~= %.12g' % (Cscaled,float(Cscaled),transfer,float(transfer)))
for fn in (SOS,RED):
    print(os.path.basename(fn),'content sha256',hashlib.sha256(read_data(fn)).hexdigest())
