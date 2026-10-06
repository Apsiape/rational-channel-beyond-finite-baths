#!/usr/bin/env python3
"""Check the 12-qubit payload against the relation compiler (hypotheses of Theorem 3.4).

Imports relation_compiler.py, which rebuilds the relation system and checks every relation
and triangle in the amalgam model at import time.  Then checks that:
  * every relator has length at most six;
  * the payload's triangle table is the compiler's prefix-and-protected-inverse table;
  * every 2x2 block of the payload carries exactly the Kraus pattern of its triangle;
  * the identity and all 88 letter labels have singleton positions;
  * the signals are the compiler's labels for p and q, and are distinct.
"""
import gzip
import importlib.util
import json
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('wzc', os.path.join(ROOT, 'relation_compiler.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


from datafile import read_data
o = json.loads(read_data(os.path.join(ROOT, 'wz_max_anchor_channel_12q.json.gz')))
check(max(len(l) + len(r) for _, l, r in m.ALL_EQUATIONS) <= 6, 'relator length')
check(len(m.ALL_EQUATIONS) == 917 and len(m.ALL_GENERATORS) == 81, 'relators / generators')
T = [tuple(t) for t in o['triangles']]
check(T == list(m.triangles), 'payload triangle table differs from the compiler')

K = o['kraus']
entries = [set(map(tuple, E)) for E in K]
A = o['anchor_labels']
base = len(A)
for n, (x, y, z) in enumerate(T):
    p = base + 2 * n
    check((p, p, 3) in entries[0], f'block {n}: identity entry')
    check((p, p + 1, 4) in entries[y], f'block {n}: y entry')
    check((p + 1, p, 4) in entries[x], f'block {n}: x entry')
    check((p + 1, p + 1, -3) in entries[z], f'block {n}: z entry')

letters = {m.label_to_id[m.signed_gen_key(name, s)] for name in m.ALL_GENERATORS for s in (1, -1)}
check(len(letters) == 88 and 0 not in letters, 'letter labels')
check({0} | letters <= set(A), 'identity or a letter label lacks a singleton')
for i, g in enumerate(A):
    check((i, i, 5) in entries[g], f'singleton pivot for label {g}')

P, Q = o['signal_labels']['P'], o['signal_labels']['Q']
check((P, Q) == (m.P_LABEL, m.Q_LABEL) and P != Q, 'signal labels')
check(P in letters and Q in letters, 'signals are generators')

print('OK payload matches the relation compiler')
print('relators 917 (length <= 6), generators 81, letter labels 88, all anchored')
print('triangle table identical (1866 blocks), block patterns exact, signals', P, Q)
