#!/usr/bin/env python3
"""Recompute the constants of Section 3 (Theorems 3.3 and 3.4), Proposition 4.1 and
Theorem 4.2 from a = 3/5 and b = 4/5, in exact rational arithmetic.

Lengths are in units of r = sqrt(d * leakage); the proof is in the paper.
"""
from fractions import Fraction as F


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


a, b = F(3, 5), F(4, 5)
check(a * a + b * b == 1, 'orthogonal block')

# Lemma 3.2: anchors.  ||B_g - V_s||_2 <= 2r; after normalization ||B_1 - I||_2 <= 2r.
anchor = F(2)

# Theorem 3.3, one prefix step x -> z = x*y with ||B_x - v(x)|| <= c r.
A_err = anchor * a + 1            # ||U_kk - aI||_2
B_err = anchor * b + 1            # ||U_{k,k+1} - b v_y||_2
check(A_err == F(11, 5) and B_err == F(13, 5), 'block entry errors')
AB = A_err + a * B_err            # ||U_kk^* U_{k,k+1} - ab v_y||_2
check(AB == F(94, 25), '94/25')
const = 1 + 1 + AB                # C-term residual + outside product + AB
check(const == F(144, 25), '144/25')
step_const = (const / b + 1) / a  # divide by b, add the e22 residual, divide by a
check(step_const == F(41, 3), '41/3')


def step(c):
    return c / a + step_const     # = (5/3) c + 41/3


cs = [anchor]
for _ in range(5):
    cs.append(step(cs[-1]))
check(cs == [F(2), F(17), F(42), F(251, 3), F(1378, 9), F(7259, 27)], 'recursion')
for ell, c in enumerate(cs, start=1):
    check(c == F(45, 2) * (F(5, 3)) ** (ell - 1) - F(41, 2), f'closed form at length {ell}')

word = cs[5] + anchor             # relator ends at label 1, ||B_1 - I|| <= 2r
check(word == F(7313, 27), 'letter-evaluation defect')
inverse = cs[1] + anchor          # protected inverse triangle has length two
check(inverse == 19, 'inverse replacement')
total = word + 6 * inverse
check(total == F(10391, 27) and total < 385 < 512, 'positive-generator defect')

# Theorem 3.4: if r < rho0/512 with rho0 <= 2, then r <= 1/256 and the coherence margin
# 7/8 - r^2 - 2r exceeds 3/4, while g = (rho0/512)^2/d <= (1/256)^2/d < 3/(4d).
r = F(1, 256)
check(F(7, 8) - r * r - 2 * r > F(3, 4), 'coherence margin')
check(r * r < F(3, 4), 'g below the coherence margin (times d)')
check(1 - F(1, 2) ** 2 / 2 == F(7, 8), 'Re tau(v_p v_q^*) >= 7/8 when ||v_p - v_q|| <= 1/2')

# Proposition 4.1: rho0 = delta/2^10.  Rounding: rho0 + 6 rho0 = 7 delta/1024 < delta/128.
check(7 * 128 < 1024, 'rounding')
# Signals: 5*2^-20 + 4 rho0 < 2^-17 since rho0 <= 2^-30 (J >= 20).
check(5 * 2 ** 10 + 4 < 2 ** 13, 'signal bound')
check(F(1, 2 ** 17) <= F(1, 2), 'beta <= 1/2')

# Theorem 4.2: g = (2^-(J+10) / 2^9)^2 / 2^12 = 2^-(2J+50).
check(2 * (10 + 9) + 12 == 50, 'gap exponent')

print('OK compiler constants')
print('94/25, 144/25, 41/3; prefix constants', ', '.join(map(str, cs)))
print('letter defect', word, '; inverse replacement', inverse, '; total', total, '< 385')
print('coherence margin > 3/(4d); g* = 2^-(2J+50)')
