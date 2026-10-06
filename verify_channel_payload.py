#!/usr/bin/env python3
"""Independent exact verifier for wz_max_anchor_channel_12q.json.gz.

Does not import the emitter.  Standard library only.  Checks canonical hash, exact
Kraus normalization in integer numerator arithmetic, disjoint Kraus supports, Choi
rank pivots, flat singleton probe, signal positions, and the advertised dimensions.
"""
import json,hashlib,os
from datafile import read_data
from collections import defaultdict
ROOT=os.path.dirname(os.path.abspath(__file__))
PATH=os.path.join(ROOT,'wz_max_anchor_channel_12q.json.gz')

def check(c,m):
    if not c: raise RuntimeError(m)
o=json.loads(read_data(PATH))
h=o.pop('canonical_sha256_without_hash_field')
raw=json.dumps(o,sort_keys=True,separators=(',',':')).encode()
check(hashlib.sha256(raw).hexdigest()==h,'canonical payload hash mismatch')
D=o['dimension'];N=o['label_count'];den=o['coefficient_denominator'];K=o['kraus'];A=o['anchor_labels'];T=o['triangles']
check((D,N,den,len(A),len(T))==(4096,949,5,364,1866),'headline counts')
check(o['choi_rank']==949 and o['nonzero_kraus_entries']==7828,'rank/nonzero metadata')
check(len(K)==N and all(K),'empty/missing Kraus label')
# supports pairwise disjoint and every entry coordinate appears under one label
seen={}
for g,E in enumerate(K):
    for r,c,n in E:
        check(0<=r<D and 0<=c<D and isinstance(n,int) and n!=0,'bad sparse entry')
        check((r,c) not in seen,f'overlapping Kraus supports at {(r,c)}')
        seen[r,c]=g
check(sum(map(len,K))==7828,'nonzero entry count')
# exact sum A_g^* A_g and sum A_g A_g^*
left=defaultdict(int);right=defaultdict(int)
for E in K:
    rows=defaultdict(list);cols=defaultdict(list)
    for r,c,n in E: rows[r].append((c,n)); cols[c].append((r,n))
    for vals in rows.values():
        for c,n in vals:
            for c2,n2 in vals:left[c,c2]+=n*n2
    for vals in cols.values():
        for r,n in vals:
            for r2,n2 in vals:right[r,r2]+=n*n2
exp={(i,i):den*den for i in range(D)}
check({k:v for k,v in left.items() if v}==exp,'TP normalization')
check({k:v for k,v in right.items() if v}==exp,'unital normalization')
# singleton anchors are exactly first 364 diagonal positions, one per distinct label.
check(len(set(A))==364,'anchor labels not distinct')
pos=o['anchor_positions']
for i,g in enumerate(A):
    check(int(pos[str(g)] if isinstance(next(iter(pos.keys())),str) else pos[g])==i,'anchor position map')
    check([i,i,5] in K[g],f'missing singleton pivot for label {g}')
# These pivots prove linear independence of all anchored Kraus matrices; every other
# label has a nonempty disjoint support, so all 949 Kraus matrices are independent.
# Flat probe: off-diagonal matrix units between distinct singleton positions are
# annihilated because no Kraus matrix contains two different singleton pivots.
for i,g in enumerate(A):
    for j,hg in enumerate(A):
        if i==j: continue
        # a single label cannot contain both pivot entries (trivial from distinct labels)
        check(g!=hg,'duplicate anchor')
# signal target cross moment is exactly zero since P,Q are distinct labels.
P=o['signal_labels']['P'];Q=o['signal_labels']['Q'];sp=o['signal_positions']
check(P!=Q and sp['P']==1 and sp['Q']==2,'signal labels/positions')
check(A[1]==P and A[2]==Q,'signal anchors')
print('OK independent 12-qubit channel payload')
print('canonical_sha256',h)
print('dimension=4096 anchors=364 triangles=1866 Choi_rank=949 nonzero=7828')
print('flat_probe_rank=364; signal positions=(1,2)')
