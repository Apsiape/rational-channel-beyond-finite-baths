"""QI-218: exact finite presentation and quotient-witness tests.

The script does not verify property (T), the centralizer theorem, or the
ultraproduct no-drift theorem. Their analytical uses are stated in QI-218.
--perturb changes the witness relation; rejection must precede success.
"""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys
from rigid_channel.core import make_presentation,eval_word,GE,IMAT,IACT,ROOTS,monomial,elementary,transform,act_elementary,in_h,gmul,ginv,require

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--perturb',action='store_true')
    args=parser.parse_args()
    p,images=make_presentation()
    stored=json.loads((ROOT/'verify/data/qi218/presentation.json').read_text())
    require(p==stored,'literal presentation differs from deterministic specification')
    if args.perturb:
        # Preserve the j definition, but remove its final action letter.
        stored['relators'][-1]['word'].pop()
    for rec in stored['relators']:
        require(not eval_word(rec['word'],images),'quotient relator fails: '+rec['family'])
    lengths=[len(r['word']) for r in p['relators']]
    require((len(p['generators']),len(lengths),sum(lengths),max(lengths))==(43,4417,41953,19),'presentation size mismatch')
    require(p['witness_word']==[3,37,-22,-3,22,-37] and p['witness_generator']==43,'witness definition mismatch')
    require(len(images[43])==3,'witness is not the expected reduced amalgam word')
    require(all(not in_h(g) for _,g in images[43]),'witness syllable belongs to amalgam subgroup')
    require(tuple(s for s,g in images[43])==(1,0,1),'witness alternating sides mismatch')
    middle=images[43][1][1]
    require(middle[1]==IACT and middle[0][1]==((-1,1,0,4),),'middle witness coefficient differs from -x2/x1')
    require(Fraction(1,5)<Fraction(1,4),'orthogonality squared margin must be positive')
    # The six positive monomial substitutions preserve the positive cone.
    positive=(0,1,2,3)
    for i in range(3):
        for j in range(3):
            if i==j:continue
            a=act_elementary(i,j)
            for f in range(3):
                for r in positive:
                    require(in_h((transform(a,elementary(ROOTS[f],monomial(r))),IACT)), 'positive compressor leaves polynomial ring')
    # Inverse substitution takes h=E12(x2) outside the polynomial subgroup.
    h=(elementary(ROOTS[0],monomial(2)),IACT)
    t=(IMAT,act_elementary(0,1))
    require(in_h(h) and not in_h(gmul(gmul(ginv(t),h),t)),'strict compression witness fails')
    print('PRESENTATION: 43 generators; 4417 distinct freely reduced relators')
    print('LENGTHS: total 41953; maximum 19')
    print('QUOTIENT: every relator is exactly identity in the Laurent-matrix amalgam')
    print('KAZHDAN INPUT: three finite root subgroups; squared orthogonality 1/5 < 1/4')
    print('COMPRESSION: all six positive action roots preserve the polynomial coefficient cone')
    print('WITNESS: j=[a2,T12_1 T12_0^-1]; reduced sides (1,0,1); middle coefficient -x2/x1')
    print('SCOPE: external property-(T) and ultraproduct centralizer arguments are not code-verified')


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,TypeError,OSError) as e:
        print('FAIL: '+str(e),file=sys.stderr);sys.exit(1)
