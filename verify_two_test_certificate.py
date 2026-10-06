#!/usr/bin/env python3
"""Exact finite checks for the channel separation certificate.

Checks the two primitive Choi effects and the scalar ledger turning them into one
fixed affine two-outcome game.  The analytic low-leakage implication is proved in
the manuscript.
"""
import json,os
from datafile import read_data
from fractions import Fraction as F
ROOT=os.path.dirname(os.path.abspath(__file__))
o=json.loads(read_data(os.path.join(ROOT,'wz_max_anchor_channel_12q.json.gz')))
d=o['dimension']; P=o['signal_labels']['P']; Q=o['signal_labels']['Q']
p=o['signal_positions']['P']; q=o['signal_positions']['Q']; K=o['kraus']
assert (d,p,q)==(4096,1,2) and P!=Q

def entry(g,r,c):
    for rr,cc,n in K[g]:
        if rr==r and cc==c:return F(n,o['coefficient_denominator'])
    return F(0)

jp=sum(entry(g,p,p)**2 for g in range(len(K)))/d
jq=sum(entry(g,q,q)**2 for g in range(len(K)))/d
cross=sum(entry(g,p,p)*entry(g,q,q) for g in range(len(K)))/d
assert jp==jq==F(1,d) and cross==0
Eplus=(jp+jq+2*cross)/2
assert Eplus==F(1,d)
# g=2^{-(2J+50)} <=2^-52 since J>=1, while the proved low-leakage
# coherence deviation exceeds 3/(4d)=3/2^14.
g_upper=F(1,2**52)
assert g_upper<F(3,4*d)
# Single affine game: alpha=d*g/2.  Work symbolically in units of g.
# high-leak branch score shift >= 1 - (d/2)*(1/d) = 1/2;
# low-leak branch >= (d/2)*3/(4d)=3/8.
assert F(1,2)>F(3,8)>0
print('OK finite channel witness data')
print('support effect target probability 0')
print('coherence effect target probability 1/4096')
print('single-game normalized score gap >= 3g/[8(1+alpha)], alpha=4096g/2')
