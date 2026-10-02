#!/usr/bin/env python3
"""Exact algebraic certificates for additional triangle pruning.

Uses no numerical optimizer. Checks implications of deleted triangle relations
by free reduction. External nonhyperlinearity/ultraproduct theorems are NOT
verified here. Run: python build_channels.py [--source PATH] [--out DIR].
"""
from __future__ import annotations
import argparse, gzip, hashlib, json
from collections import Counter
from fractions import Fraction
from pathlib import Path

EXPECTED = '38150decf3d203e3eabbc5dee8d8a3f2d6390aa960c50f90fe4c87b804729b68'

def need(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)

def red(word):
    result = []
    for x in word:
        if result and result[-1] == -x:
            result.pop()
        else:
            result.append(x)
    return result

def inverse(word):
    return [-x for x in reversed(word)]

def symbol(g):
    # Label zero is the identity, not a free generator.
    return [] if g == 0 else [g + 1]

def word(t):
    x, y, z = t
    return red(symbol(x) + symbol(y) + inverse(symbol(z)))

def main(source: Path, out: Path) -> None:
    raw = gzip.decompress(source.read_bytes())
    need(hashlib.sha256(raw).hexdigest() == EXPECTED, 'unexpected baseline fingerprint')
    c = json.loads(raw)
    T = list(map(tuple, c['triangles']))
    lookup = {t: i for i, t in enumerate(T)}
    need(len(lookup) == len(T), 'baseline triangles not unique')
    G = dict(c['generator_labels'])
    inv = {G[s]: G[-s] for s in G}
    foundation = {(G[s], G[-s], 0) for s in range(1, 44)}
    need(foundation <= set(T), 'missing positive-generator inverse equations')
    inv_reference = {}
    for s in range(1, 44):
        t = (G[s], G[-s], 0)
        inv_reference[G[s]] = (lookup[t], False)
        inv_reference[G[-s]] = (lookup[t], True)
    certificates = []
    deleted = set()
    for t in T:
        x, y, z = t
        if y not in inv:
            continue
        r = (z, inv[y], x)
        if r not in lookup or not t < r:
            continue
        # Keep t, delete r; foundation triangles cannot occur in these pairs.
        need(t not in foundation and r not in foundation, 'unexpected foundation overlap')
        f, reversed_orientation = inv_reference[y]
        conjugator = symbol(x) + (symbol(y) if reversed_orientation else [])
        certificates.append({'deleted': lookup[r], 'kind': 'right_inverse',
            'cells': [{'reference': lookup[t], 'sign': -1, 'conjugator': []},
                      {'reference': f, 'sign': 1, 'conjugator': conjugator}]})
        deleted.add(lookup[r])
    right_count = len(deleted)
    for t in T:
        x, y, z = t
        if z != 0 or t in foundation:
            continue
        r = (y, x, 0)
        need(r in foundation, 'unhandled zero-output relation')
        need(lookup[t] not in deleted, 'duplicate deletion')
        certificates.append({'deleted': lookup[t], 'kind': 'reversed_inverse',
            'cells': [{'reference': lookup[r], 'sign': 1, 'conjugator': symbol(x)}]})
        deleted.add(lookup[t])
    retained = [i for i in range(len(T)) if i not in deleted]
    keep = set(retained)
    for cert in certificates:
        product = []
        for cell in cert['cells']:
            need(cell['reference'] in keep, 'certificate depends on deleted equation')
            r = word(T[cell['reference']])
            if cell['sign'] == -1:
                r = inverse(r)
            h = cell['conjugator']
            product.extend(h + r + inverse(h))
        need(red(product) == word(T[cert['deleted']]), 'free-group certificate failed')
    all_labels = set(range(len(c['labels'])))
    covered = {g for i in retained for g in T[i]}
    need(all_labels <= covered, 'unitary-forcing coverage lost')
    out.mkdir(parents=True, exist_ok=True)
    cert_data = {'baseline_sha256': EXPECTED, 'original_triangle_count': len(T),
                'retained_original_indices': retained,
                'free_generator_convention': 'label g>0 is free generator g+1; label 0 is identity',
                'certificates': certificates}
    (out/'deletion_certificates.json').write_text(json.dumps(cert_data, separators=(',', ':'))+'\n', newline='\n')
    summary = {'right_inverse_pair_deletions': right_count,
               'reversed_generator_inverse_deletions': len(deleted)-right_count,
               'total_deleted': len(deleted), 'retained_triangles': len(retained),
               'certificate_count': len(certificates), 'all_labels_covered': True,
               'variants': {},
               'scope': 'Exact finite data and free-group implications only. Global nonmembership inherits baseline analytical inputs; global numerical gap remains unevaluated.'}
    for name, anchors, padded in [('one_anchor', [0], None),
                                  ('two_anchor', [0, G[43]], None),
                                  ('fifteen_qubits', [0, G[43]], 32768)]:
        natural_d = len(anchors) + 2*len(retained)
        d = natural_d if padded is None else padded
        need(d >= natural_d, 'padding dimension too small')
        kraus = [[] for _ in c['labels']]
        for p, g in enumerate(anchors):
            kraus[g].append([p, p, 5])
        for n, index in enumerate(retained):
            x,y,z = T[index]; p = len(anchors)+2*n
            for g,entry in [(0,[p,p,3]),(y,[p,p+1,4]),(x,[p+1,p,4]),(z,[p+1,p+1,-3])]:
                kraus[g].append(entry)
        for p in range(natural_d, d):
            kraus[0].append([p,p,5])
        rows=[0]*d; cols=[0]*d; slots=set(); weights=[]
        for entries in kraus:
            need(bool(entries), 'empty label')
            gr=set();gc=set();weight=0
            for a,b,n in entries:
                need((a,b) not in slots, 'Kraus supports overlap')
                need(a not in gr and b not in gc, 'Kraus not a weighted partial permutation')
                slots.add((a,b));gr.add(a);gc.add(b)
                rows[a]+=n*n;cols[b]+=n*n;weight+=n*n
            weights.append(weight)
        need(rows==[25]*d and cols==[25]*d, 'normalization failure')
        obj={'format':'rational-triangle-certified-pruning-v1',
             'entry_denominator':5, 'dimension':d, 'kraus_rank':len(kraus),
             'anchor_labels':anchors, 'identity_padding_count':d-natural_d,
             'labels':c['labels'], 'generator_labels':c['generator_labels'],
             'triangles':[list(T[i]) for i in retained], 'kraus':kraus,
             'baseline_uncompressed_sha256':EXPECTED}
        data=json.dumps(obj,separators=(',',':'),sort_keys=True).encode()
        (out/f'{name}_channel.json.gz').write_bytes(gzip.compress(data,mtime=0))
        summary['variants'][name]={'dimension':d,'rank':len(kraus),
             'sparse_entries':len(slots),'identity_padding_count':d-natural_d,
             'identity_weight_numerator':weights[0], 'j_weight_numerator':weights[G[43]],
             'restricted_choi_gap_squared':str(Fraction(weights[0]*weights[G[43]],(25*d)**2)),
             'global_numerical_gap':None, 'sha256_uncompressed':hashlib.sha256(data).hexdigest()}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n', newline='\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    here=Path(__file__).resolve().parent
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=here/'baseline/verify/data/qi219/channel.json.gz')
    p.add_argument('--out',type=Path,default=here/'data')
    args=p.parse_args()
    main(args.source,args.out)
