#!/usr/bin/env python3
"""Recompute the constants of Theorems 7.1 and 7.2, Lemma 7.3 and the Section 5 witness margin.

Lengths are in units of r = sqrt(d * leakage) <= 1; the proof is in the paper.
"""
from fractions import Fraction as F


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


a, b = F(3, 5), F(4, 5)

# Letters from triangles: ||U*U - a^2||_2 <= (1 + a)(2a + 1) r = 88/25 r; outside part <= r.
uu = (1 + a) * (2 * a + 1)
check(uu == F(88, 25), '88/25')
yy = (uu + 1) / (b * b)           # ||Y*Y - I||_2, and ||Y - V||_2 <= ||Y*Y - I||_2
check(yy == F(113, 16), '113/16')
letter = yy + 1 / b               # plus ||B_x - Y||_2 <= r/b
check(letter == F(133, 16), 'letters within 133/16 r')


def solve(c):
    """The bound in the proof of Theorem 7.1 on ||B_z - v_x v_y|| / r when x, y are within c r."""
    return ((3 + 3 * a + b * c * (1 + a)) / b + 1) / a


z_only = solve(letter)
check(z_only == F(203, 6) and z_only < 34, 'z-only labels within 203/6 r < 34 r')
ret = solve(z_only) + z_only      # conservative: x, y within the z-only bound
check(ret == F(2443, 18) and ret < 136, 'retained triangles below 136 r')
check(2 * 136 == 272, 'deleted triangles below 272 r')
rel = 18 * 272 + 19 * 136
check(rel == 7480 and rel < 2 ** 13, 'relator ledger below 2^13 r')

# Theorem 7.2.  g = delta^2 / (2^26 d); nu < g gives r^2 = d nu < delta^2 / 2^26, so r < 2^-13 delta <= 2^-13.
check(F(26, 2) == 13, 'r < 2^-13 delta')
check(2 * 13 + 15 == 41, 'gap exponent 2J + 41 at d <= 2^15')
check(F(1, 2 ** 13) * 34 <= F(1, 2), '||B_j|| <= 1 + 34 r <= 3/2')
check(2 * F(3, 2) + 34 == 37, '||B_1 - I|| ||B_j|| + ||B_j - v_j|| <= 2r * 3/2 + 34 r = 37 r')
check(1 - F(1, 2 ** 17) ** 2 / 2 == 1 - F(1, 2 ** 35), 'Re tau(v_j) >= 1 - 2^-35')
check(a * a <= F(16, 25), 'sqrt(s_1 s_j) >= 3/5 needs s_j >= 9/25')
d_max, r_max = 2 ** 15, F(1, 2 ** 60)
check(F(3, 5) * (1 - F(1, 2 ** 35) - 37 * r_max) - 34 * r_max * d_max > F(1, 2), 'coherence test exceeds 1/(2d)')

# Witness for the reduced channels (Section 5): s_1 + s_j < 2^13 from the shipped Kraus data, so the
# margin g / (2 (s_1 + s_j)(1 + alpha)) exceeds g / 2^15 (alpha <= 1).
import json, os
from datafile import read_data
HERE = os.path.dirname(os.path.abspath(__file__))
for name in ['fifteen_qubits_channel', 'two_anchor_channel', 'one_anchor_channel']:
    o = json.loads(read_data(os.path.join(HERE, 'earlier', name + '.json.gz')))
    weight = lambda g: F(sum(v * v for _, _, v in o['kraus'][g]), 25)
    s1, sj = weight(0), weight(85)
    check(s1 >= 1 and sj >= F(9, 25), f'{name}: singleton weights')
    check(sum(weight(g) for g in range(len(o['kraus']))) == o['dimension'], f'{name}: sum of s_g = d')
    check(s1 + sj < 2 ** 13 and 2 * (s1 + sj) * 2 < 2 ** 15, f'{name}: witness margin above g / 2^15')
    print(f'{name}: s_1 + s_j = {s1 + sj}')

# Lemma 7.3 (earlier report, Theorem 6).  Entries: sqrt(1+c^2) + |c| < 3, and sqrt(2(1+c^2)) + |c| < 3.
for c in (a, b):
    check(1 + c * c < (3 - c) ** 2 and 2 * (1 + c * c) < (3 - c) ** 2, 'entry constant 3k')
check(13 * F(25, 12) < 32 and (12 + F(3, 2)) * F(25, 12) < 32, 'triangle defect below 32 k (sqrt 2 < 3/2)')
check(64 * 2 == 128, 'k <= e_0 / (64 l_*) <= 1/128')
check(F(1, 128) ** 2 * 2 + 2 * F(1, 128) < F(7, 8), 'witness trace contradiction')
check(4096 * 19 ** 2 == 1478656 and 1478656 * 44086 == 65188028416, 'Corollary 7.4 denominators')
check(2 ** 35 < 65188028416 < 2 ** 36 and 2 ** 20 < 1478656 < 2 ** 21, 'full-channel exponents')

print('OK earlier-channel leakage constants')
print('letters within', letter, 'r; z-only within', z_only, 'r; retained', ret, '< 136 r')
print('deleted < 272 r; relators < 7480 r < 2^13 r; gap 2^-(2J+41)')
print('fully referenced channel: Choi > 2^-(2J+36), diamond > 2^-(2J+21)')
