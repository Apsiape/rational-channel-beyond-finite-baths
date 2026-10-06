# Adapted from the verification package of the earlier report (Mghirbi and Douglas, Zenodo,
# doi:10.5281/zenodo.23093839; MIT licence, see LICENSE-MIT.txt), with paths set to this release.
"""Independently evaluate all labels and multiplication triangles of the fully referenced channel.

The complete sparse rational channel is read from channel.json.gz.  Checks
reconstruct the generator images, relator paths and rational channel identities;
none certifies the external matrix-ultraproduct theorem by finite sampling.
--perturb corrupts a Kraus coefficient. --perturb-label corrupts a label moment.
"""
import argparse
from collections import Counter,defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import sys
from rigid_channel.independent import (need,generator_images,eval_word,decode,equal,multiply,normal,hfactor)

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--perturb',action='store_true')
    ap.add_argument('--perturb-label',action='store_true')
    args=ap.parse_args()
    p=json.loads((ROOT/'presentation.json').read_text())
    gz=ROOT/'original_channel.json.gz'
    raw=gzip.decompress(gz.read_bytes()) if gz.exists() else (ROOT/'original_channel.json').read_bytes()
    c=json.loads(raw)
    need(len(p['generators'])==43 and len(p['relators'])==4417,'literal presentation size')
    need(p['witness_generator']==43 and p['witness_word']==[3,37,-22,-3,22,-37],'distinguished word binding')
    need(p['relators'][-1]['word']==[-43,3,37,-22,-3,22,-37],'distinguished relator binding')
    if args.perturb:c['kraus'][0][0][2]=4
    if args.perturb_label:c['labels'][1]['word']=[2]
    gens=generator_images(p)
    for r in p['relators']:
        need(not eval_word(r['word'],gens),'independent relator evaluation failed')
    vals=[];simple={};multi=[]
    for idx,label in enumerate(c['labels']):
        need(label['id']==idx,'label order mismatch')
        value=decode(label['normal_form'])
        need(equal(value,eval_word(label['word'],gens)),'label is not its declared group word')
        need(value==normal(value),'label is not reduced')
        if len(value)<=1:
            need(value not in simple,'duplicate single-factor label');simple[value]=idx
        else:
            need(all(not hfactor(f) for f in value),'reduced word contains polynomial factor')
            need(all(value[k][0]!=value[k+1][0] for k in range(len(value)-1)),'nonalternating label')
            for k in multi:need(not equal(value,vals[k]),'duplicate multi-factor label')
            multi.append(idx)
        vals.append(value)
    need(not vals[0],'anchor zero is not identity')
    L=len(vals);T=len(c['triangles']);d=c['dimension']
    need((d,L,T)==(44086,11230,16428),'channel dimensions disagree')
    gl=dict(c['generator_labels'])
    need(set(gl)==set(gens),'generator label domain disagrees')
    for k,g in gens.items():need(equal(g,vals[gl[k]]),'wrong generator label')
    need(gl[43]!=0 and gl[-43]!=0,'nontrivial witness anchor was lost')
    table={}
    for x,y,z in c['triangles']:
        need(x!=0 and y!=0,'trivial triangle retained')
        need((x,y) not in table,'duplicate multiplication row')
        table[x,y]=z
        need(equal(multiply(vals[x],vals[y]),vals[z]),'false multiplication triangle')
    def step(x,y):
        if x==0:return y
        if y==0:return x
        need((x,y) in table,'missing multiplication relation')
        return table[x,y]
    for r in p['relators']:
        acc=0
        for k in r['word']:acc=step(acc,gl[k])
        need(acc==0,'relator does not close in channel triangles')
    for k in range(1,44):need(step(gl[k],gl[-k])==0,'missing inverse triangle')
    expected=[[[k,k,5]] for k in range(L)]
    for n,(x,y,z) in enumerate(c['triangles']):
        i=L+2*n
        expected[0].append([i,i,3]);expected[y].append([i,i+1,4])
        expected[x].append([i+1,i,4]);expected[z].append([i+1,i+1,-3])
    need(c['entry_denominator']==5 and c['kraus_rank']==L,'channel encoding changed')
    need(c['kraus']==expected,'Kraus data are not the anchored rational triangle construction')
    row_norm=[0]*d;col_norm=[0]*d;slots=set();weights=[]
    for k,entries in enumerate(c['kraus']):
        rows=set();cols=set();weight=0
        for i,j,n in entries:
            need(0<=i<d and 0<=j<d and n in (5,3,4,-3),'invalid sparse entry')
            need((i,j) not in slots,'overlap between Kraus supports');slots.add((i,j))
            need(i not in rows and j not in cols,'Kraus matrix is not a weighted partial permutation')
            rows.add(i);cols.add(j);row_norm[i]+=n*n;col_norm[j]+=n*n;weight+=n*n
        need(entries[0]==[k,k,5],'missing independent singleton anchor')
        weights.append(weight)
    need(row_norm==[25]*d and col_norm==[25]*d,'unital or trace-preserving identity fails')
    need(len(slots)==76942 and sum(weights)==25*d,'slot/spectral normalization mismatch')
    spectrum=Counter(weights)
    info=json.loads((ROOT/'model/instance.json').read_text())
    need(info['presentation_sha256']==hashlib.sha256((ROOT/'presentation.json').read_bytes()).hexdigest(),'presentation fingerprint mismatch')
    need(info['channel_uncompressed_sha256']==hashlib.sha256(raw).hexdigest(),'channel byte fingerprint mismatch')
    need(info['choi_weight_histogram']==[[k,spectrum[k]] for k in sorted(spectrum)],'Choi spectral weight histogram changed')
    need(info['numerical_gap_lower_bound'] is None,'unevaluated distance gap must not be a numerical certificate')
    print('INDEPENDENT: 4417 relators, 11230 distinct labels and 16428 triangles checked exactly')
    print('ENCODING: every defining word and generator inverse closes using retained triangles')
    print('CHANNEL: dimension 44086; Choi rank 11230; 76942 sparse rational Kraus entries')
    print('COEFFICIENTS: denominator 5; numerators 5, 3, 4, -3')
    print('CPTP: sum K_g* K_g = I; UNITAL: sum K_g K_g* = I, by integer identities')
    print('CHOI: disjoint supports; positive weights; exact trace one; singleton rank-11230 flat probe')
    print('SCOPE: labels, triangles and channel identities only; the distance bounds are in the paper')


if __name__=='__main__':
    try:main()
    except (ValueError,KeyError,TypeError,OSError) as e:
        print('FAIL: '+str(e),file=sys.stderr);sys.exit(1)
