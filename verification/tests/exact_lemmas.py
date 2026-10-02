#!/usr/bin/env python3
"""Exact arithmetic behind numerical constants and the two-position probe.

The universal analytical estimates themselves are proved in the manuscripts.
"""
from __future__ import annotations
from fractions import Fraction as F
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from channel_api import RationalChannel


def need(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def main() -> None:
    checked = []
    for c in (F(3, 5), F(4, 5)):
        need(2 * (1 + c*c) < (3 - c)**2, 'two-input extraction coefficient')
    need((12 + F(3, 2)) * F(25, 12) < 32 and F(3, 2)**2 > 2,
         'row-product extraction coefficient')
    need(13 * F(25, 12) < 32, 'Choi extraction coefficient')
    need(4096 * 19**2 == 1478656, 'half-diamond denominator')
    need(4096 * 19**2 * 44086 == 65188028416, 'Choi denominator')
    need(2 - 2*(2*F(1,4)**2 + 2*F(1,4)) == F(3,4) > F(1,4),
         'trace contradiction')
    need(F(112,24) == F(14,3), 'memory leading coefficient ratio')
    # ln(2) > 9/13 follows from the first three positive terms of
    # 2 * sum_{k>=0} (1/3)^(2k+1)/(2k+1).
    lower_ln2 = 2 * (F(1,3) + F(1,81) + F(1,1215))
    need(lower_ln2 > F(9,13), 'purity-budget logarithm bound')
    checked += ['conditional extraction constants', 'memory coefficient', 'purity-budget logarithm bound']
    for m in (2, 10, 100, 1000):
        real, imag = F(m*m-1,m*m+1), F(2*m,m*m+1)
        need(real*real+imag*imag == 1, 'rational phase unitarity')
        need((1-real)**2+imag**2 == F(4,m*m+1), 'phase displacement')
        need(((1-real)**2+imag**2)/172 == F(1,43*(m*m+1)), 'conjugation energy')
    checked.append('four exact rational-phase examples')
    for variant in ('two_anchor','fifteen_qubits'):
        channel = RationalChannel(ROOT/'data'/f'{variant}_channel.json.gz')
        need(channel.image_of_matrix_unit(0,0) == {(0,0):F(1)}, 'first singleton')
        need(channel.image_of_matrix_unit(1,1) == {(1,1):F(1)}, 'witness singleton')
        need(channel.image_of_matrix_unit(0,1) == {}, 'dephasing cross term')
        plus = [(i,j,F(1,2)) for i in (0,1) for j in (0,1)]
        need(channel.apply_sparse(plus) == {(0,0):F(1,2),(1,1):F(1,2)}, 'flat rank-two probe')
        # The difference between |+><+| and I_2/2 squares to I_2/4.
        difference = [[F(0),F(1,2)],[F(1,2),F(0)]]
        square = [[sum(difference[i][k]*difference[k][j] for k in (0,1))
                   for j in (0,1)] for i in (0,1)]
        need(square == [[F(1,4),F(0)],[F(0),F(1,4)]], 'half-distance calculation')
    checked.append('two-position probes and half-distance 1/2')
    print(json.dumps({'passed': True, 'checked': checked,
                      'scope': 'Exact scalar and sparse identities, not a proof of general rigidity.'}, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, IndexError, OSError) as exc:
        print('FAIL: '+str(exc), file=sys.stderr)
        sys.exit(1)
