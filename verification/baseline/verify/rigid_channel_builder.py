"""QI-219: rebuild the actual rational channel, without an upstream compiler.

Default: regenerate and compare with the retained exact presentation/channel.
--emit DIR writes presentation.json, channel.json.gz and instance.json to DIR.
Neither mode computes a numerical lower bound on the tracial distance.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
from rigid_channel.core import make_presentation,compile_gadget,require

ROOT=Path(__file__).resolve().parents[1]


def canonical(obj):
    return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()


def build():
    p,images=make_presentation();c=compile_gadget(p,images)
    raw_p=canonical(p);raw_c=canonical(c)
    weights=Counter(sum(v*v for _,_,v in k) for k in c['kraus'])
    summary={
        'format':'qi219-instance-summary-v1','channel_file':'channel.json.gz',
        'channel_uncompressed_sha256':hashlib.sha256(raw_c).hexdigest(),
        'presentation_sha256':hashlib.sha256(raw_p).hexdigest(),
        'system_dimension':c['dimension'],'choi_rank':c['kraus_rank'],
        'sparse_kraus_entries':sum(map(len,c['kraus'])),
        'choi_weight_denominator':25*c['dimension'],
        'choi_weight_histogram':[[a,weights[a]] for a in sorted(weights)],
        'numerical_gap_lower_bound':None,
        'gap_status':'Positive by the QI-218/219 analytical argument; not evaluated numerically.',
        'external_inputs':['arXiv:0809.4095v2 Theorems 1.2 and 6.2',
                           'arXiv:2609.37761v1 Theorem B',
                           'arXiv:2609.27795v1 Theorem 1.2'],
        'claim_status':'CANDIDATE'}
    return p,c,summary,raw_p,raw_c


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--emit',type=Path,metavar='DIRECTORY')
    ap.add_argument('--perturb',action='store_true')
    args=ap.parse_args()
    p,c,s,rp,rc=build()
    if args.perturb:
        c['triangles'][0][2]=0
        rc=canonical(c)
    if args.emit is not None:
        require(not args.perturb,'do not emit deliberately corrupted data')
        args.emit.mkdir(parents=True,exist_ok=True)
        (args.emit/'presentation.json').write_bytes(rp)
        (args.emit/'channel.json.gz').write_bytes(gzip.compress(rc,mtime=0))
        (args.emit/'instance.json').write_text(json.dumps(s,indent=2)+'\n', newline='\n')
    else:
        require(rp==(ROOT/'verify/data/qi218/presentation.json').read_bytes(),'recompiled presentation bytes differ')
        stored=gzip.decompress((ROOT/'verify/data/qi219/channel.json.gz').read_bytes())
        require(rc==stored,'recompiled channel bytes differ')
        require(s==json.loads((ROOT/'verify/data/qi219/instance.json').read_text()),'instance summary differs')
    require((len(p['generators']),len(p['relators']))==(43,4417),'presentation dimensions')
    require((c['dimension'],c['kraus_rank'],len(c['triangles']))==(44086,11230,16428),'channel dimensions')
    require(sum(p['raw_relation_counts'].values())==4831,'raw relator total')
    require(4831-len(p['relators'])==414,'common positive relator total')
    require(s['numerical_gap_lower_bound'] is None,'distance has not been evaluated')
    print('BUILD: complete finite presentation and rational channel regenerated')
    print('COUNTS: 43 generators; 4417 relators; 11230 labels; 16428 triangles; dimension 44086')
    print('PRESENTATION SHA256: '+s['presentation_sha256'])
    print('CHANNEL JSON SHA256: '+s['channel_uncompressed_sha256'])
    print('SOURCE OF COEFFICIENTS: explicit finite Laurent-matrix arithmetic, not a MIP* compiler')
    print('GAP: strictly positive by the analytical argument; numerical lower bound remains null')


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,TypeError,OSError) as e:
        print('FAIL: '+str(e),file=sys.stderr);sys.exit(1)
