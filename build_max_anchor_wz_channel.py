#!/usr/bin/env python3
"""Build/verify the 12-qubit maximal-distinct-anchor WZ relation channel.

The 1866 triangle blocks consume 3732 system coordinates. Exactly 364 singleton
coordinates remain in dimension 4096. We assign them to 364 distinct labels,
including 0,P,Q, maximizing the immediately visible flat dephasing probe among
12-qubit channels with this fixed triangle table.
"""
import importlib.util, pathlib, gzip, json, hashlib
from collections import defaultdict
HERE=pathlib.Path(__file__).resolve().parent
SRC=HERE/'relation_compiler.py'
spec=importlib.util.spec_from_file_location('wzc',SRC)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
triangles=m.triangles; N=m.N_LABELS; P=m.P_LABEL; Q=m.Q_LABEL
assert (N,len(triangles))==(949,1866)
D=4096
anchor_count=D-2*len(triangles)
assert anchor_count==364
anchors=[0,P,Q]+[g for g in range(N) if g not in {0,P,Q}][:anchor_count-3]
assert len(anchors)==len(set(anchors))==364
covered=set(anchors)
for x,y,z in triangles: covered.update((0,x,y,z))
assert covered==set(range(N)), (len(covered),N)
kraus=[[] for _ in range(N)]
for pos,g in enumerate(anchors): kraus[g].append([pos,pos,5])
for n,(x,y,z) in enumerate(triangles):
    p=len(anchors)+2*n
    kraus[0].append([p,p,3]); kraus[y].append([p,p+1,4])
    kraus[x].append([p+1,p,4]); kraus[z].append([p+1,p+1,-3])
assert len(anchors)+2*len(triangles)==D
# exact TP and unital normalization in numerator/5 arithmetic
left=defaultdict(int); right=defaultdict(int)
for E in kraus:
    br=defaultdict(list); bc=defaultdict(list)
    for r,c,n in E: br[r].append((c,n)); bc[c].append((r,n))
    for vals in br.values():
        for c1,n1 in vals:
            for c2,n2 in vals: left[c1,c2]+=n1*n2
    for vals in bc.values():
        for r1,n1 in vals:
            for r2,n2 in vals: right[r1,r2]+=n1*n2
expected={(i,i):25 for i in range(D)}
assert {k:v for k,v in left.items() if v}==expected
assert {k:v for k,v in right.items() if v}==expected
assert all(kraus[g] for g in range(N))
nonzero=sum(map(len,kraus)); assert nonzero==4*len(triangles)+len(anchors)==7828
positions={str(g):i for i,g in enumerate(anchors)}
# The regular-trace target is complete dephasing on distinct singleton labels.
# Thus maximally entangled input supported on these positions has a flat rank-364 output.
obj={
 'name':'wz_max_anchor_channel_12q','dimension':D,'qubits':12,
 'coefficient_denominator':5,'label_count':N,'choi_rank':N,
 'triangle_count':len(triangles),'anchor_count':len(anchors),
 'anchor_labels':anchors,'anchor_positions':positions,
 'signal_labels':{'P':P,'Q':Q},'signal_positions':{'P':positions[str(P)],'Q':positions[str(Q)]},
 'flat_probe_rank':len(anchors),'natural_dimension':D,
 'nonzero_kraus_entries':nonzero,'triangles':triangles,'kraus':kraus,
 'gap_formulas':{
   'J':'2^(2^1122000)','delta_rigidity':'2^-J',
   'compiled_relator_tolerance':'2^(-(J+10))',
   'support_extraction_r0':'2^(-(J+19))',
   'global_choi_gap':'2^(-(2J+50))',
   'global_half_diamond_gap':'2^(-(2J+50))',
   'finite_two_test_certificate_gap':'2^(-(2J+50))'
 }
}
raw=json.dumps(obj,sort_keys=True,separators=(',',':')).encode()
sha=hashlib.sha256(raw).hexdigest(); obj['canonical_sha256_without_hash_field']=sha
import sys
out=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'wz_max_anchor_channel_12q.json.gz'
with open(out,'wb') as f:
    with gzip.GzipFile(filename='',mode='wb',fileobj=f,mtime=0) as gz:
        gz.write(json.dumps(obj,sort_keys=True,separators=(',',':')).encode())
print('OK maximal-anchor channel')
print('dimension',D,'anchors',len(anchors),'triangles',len(triangles),'labels/rank',N,'nonzero',nonzero)
print('P,Q positions',positions[str(P)],positions[str(Q)])
print('flat_probe_rank',len(anchors))
print('sha256',sha)
print(out)
