#!/usr/bin/env python3
"""Emit two anchor-pruned channels from retained QI-219 finite data.

No numerical optimizer or external group compiler is required. Nonmembership
still depends on the baseline analytical theorem; this program proves no
matrix-dimension-independent global gap.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json
from fractions import Fraction
from pathlib import Path
ROOT = Path(__file__).resolve().parent
BASE = ROOT / 'verify/data/qi219/channel.json.gz'
EXPECTED = '38150decf3d203e3eabbc5dee8d8a3f2d6390aa960c50f90fe4c87b804729b68'
def emit(destination: Path) -> None:
    raw = gzip.decompress(BASE.read_bytes())
    if hashlib.sha256(raw).hexdigest() != EXPECTED:
        raise ValueError('baseline fingerprint differs; run semantic baseline verification')
    old = json.loads(raw)
    j = dict(old['generator_labels'])[43]
    destination.mkdir(parents=True, exist_ok=True)
    for name, anchors, bound in (
        ('one_anchor', [0], Fraction(1, 533)),
        ('two_anchor', [0, j], Fraction(3, 1000)),
    ):
        c = {k: v for k, v in old.items() if k not in ('dimension', 'format', 'kraus')}
        c['format'] = 'rational-triangle-pruned-anchors-v1'
        c['anchor_labels'] = anchors
        c['dimension'] = len(anchors) + 2 * len(old['triangles'])
        c['kraus'] = [[] for _ in old['labels']]
        for row, g in enumerate(anchors):
            c['kraus'][g].append([row, row, 5])
        for n, (x, y, z) in enumerate(old['triangles']):
            p = len(anchors) + 2*n
            for g, entry in ((0,[p,p,3]), (y,[p,p+1,4]),
                             (x,[p+1,p,4]), (z,[p+1,p+1,-3])):
                c['kraus'][g].append(entry)
        weights = [sum(v*v for _,_,v in entries) for entries in c['kraus']]
        squared = Fraction(weights[0]*weights[j], (25*c['dimension'])**2)
        if not bound*bound < squared:
            raise ValueError('claimed restricted bound exceeds the exact certificate')
        c['baseline_uncompressed_sha256'] = EXPECTED
        out = json.dumps(c, separators=(',',':'), sort_keys=True).encode()
        (destination/f'{name}_channel.json.gz').write_bytes(
            gzip.compress(out, compresslevel=9, mtime=0))
        summary = {
            'variant': name, 'anchor_labels': anchors,
            'dimension': c['dimension'], 'rank':len(c['kraus']),
            'sparse_entries':sum(map(len,c['kraus'])),
            'identity_weight_numerator':weights[0], 'witness_label':j,
            'witness_weight_numerator':weights[j],
            'spectrum_denominator':25*c['dimension'],
            'restricted_gap_squared':str(squared),
            'restricted_support_gap_lower_bound':str(bound),
            'global_gap_lower_bound':None,
            'restricted_comparison_class':
              'finite matrix-tracial factorizations with Choi support contained in the target Choi support',
            'nonmembership_status':
              'inherits the baseline candidate analytic argument; not certified by finite checks',
            'sha256_uncompressed':hashlib.sha256(out).hexdigest(),
        }
        (destination/f'{name}_summary.json').write_text(json.dumps(summary,indent=2)+'\n', newline='\n')
        print(f'{name}: d={c["dimension"]}, rank={len(c["kraus"])}, restricted gap >= {bound}')
    print('Global numerical gap: unresolved. No general rigidity theorem is code-certified.')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--emit', type=Path, default=ROOT/'data')
    args = parser.parse_args()
    emit(args.emit)
