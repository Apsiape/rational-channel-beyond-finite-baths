#!/usr/bin/env python3
"""Check pruning by remapping old slots, independently of the new emitter.

Baseline group and word checks are separate prerequisites. This program checks
finite numerical data and scoped certificate arithmetic, not the general proofs.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def need(test, message):
    if not test: raise ValueError(message)
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--variant', choices=['one_anchor','two_anchor'],default='two_anchor')
    ap.add_argument('--perturb',choices=['coefficient','triangle','scope'])
    args=ap.parse_args()
    old=json.loads(gzip.decompress((ROOT/'verify/data/qi219/channel.json.gz').read_bytes()))
    raw=gzip.decompress((ROOT/f'data/{args.variant}_channel.json.gz').read_bytes())
    c=json.loads(raw); info=json.loads((ROOT/f'data/{args.variant}_summary.json').read_text())
    if args.perturb=='coefficient': c['kraus'][0][0][2]=4
    if args.perturb=='triangle': c['triangles'].pop()
    if args.perturb=='scope': info['global_gap_lower_bound']='3/1000'
    j=dict(old['generator_labels'])[43]
    anchors=[0] if args.variant=='one_anchor' else [0,j]
    need(c['anchor_labels']==anchors,'wrong anchors')
    for field in ('labels','generator_labels','triangles'):
        need(c[field]==old[field],f'{field} changed during anchor pruning')
    L=len(old['labels']); m=len(old['triangles']); d=len(anchors)+2*m
    need(c['dimension']==d and len(c['kraus'])==L,'dimension or rank mismatch')
    mapping={g:i for i,g in enumerate(anchors)}
    def new_index(i):
        return mapping.get(i) if i<L else i-L+len(anchors)
    for g, entries in enumerate(old['kraus']):
        mapped=[]
        for a,b,v in entries:
            x,y=new_index(a),new_index(b)
            if x is not None and y is not None: mapped.append([x,y,v])
        need(c['kraus'][g]==mapped,'pruned Kraus entries differ from exact old-slot remapping')
    seen={z for tri in c['triangles'] for z in tri}
    need(set(range(L))<=seen|set(anchors),'a label has no unitary-forcing occurrence')
    row=[0]*d;col=[0]*d;slots=set();weights=[]
    for entries in c['kraus']:
        need(entries,'empty Kraus operator would invalidate retained rank')
        rows=set();cols=set();weight=0
        for a,b,v in entries:
            need(0<=a<d and 0<=b<d,'bad entry index')
            need((a,b) not in slots,'overlapping Kraus supports')
            need(a not in rows and b not in cols,'Kraus operator is not a partial permutation')
            slots.add((a,b));rows.add(a);cols.add(b)
            row[a]+=v*v;col[b]+=v*v;weight+=v*v
        weights.append(weight)
    need(row==[25]*d and col==[25]*d,'Kraus normalization failed')
    need(sum(weights)==25*d,'Choi trace is not one')
    squared=Fraction(weights[0]*weights[j],(25*d)**2)
    bound=Fraction(info['restricted_support_gap_lower_bound'])
    need(Fraction(info['restricted_gap_squared'])==squared,'wrong witness contraction')
    need(bound>0 and bound*bound<squared,'restricted rational bound failed')
    need(info['global_gap_lower_bound'] is None,'a support-restricted bound was mislabelled global')
    need(info['sha256_uncompressed']==hashlib.sha256(raw).hexdigest(),'binding mismatch')
    need(info['sparse_entries']==len(slots) and info['rank']==L,'wrong finite summary')
    print(f'PASS {args.variant}: d={d}; rank={L}; nonzeros={len(slots)}; exact CPTP/unital identities')
    print('PASS tables and relator paths unchanged; every label has a unitary-forcing occurrence')
    print(f'PASS support-restricted gap arithmetic: squared={squared}; rational bound={bound}')
    print('SCOPE: global numerical gap unresolved; general rigidity not formally certified')
if __name__=='__main__':
    try: main()
    except (ValueError,KeyError,TypeError,OSError) as e:
        print('FAIL: '+str(e),file=sys.stderr);sys.exit(1)
