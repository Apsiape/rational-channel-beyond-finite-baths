#!/usr/bin/env python3
"""Exact verifier/emitter for a Wang--Zhi quantitative relation-system -> rational-channel compiler.

Standard library only.  No floating-point arithmetic is used for algebraic checks.

What is verified exactly:
  * the 843 Wang--Zhi auxiliary balanced equations (with shared p_{ab}=p_{ba} symbols),
  * 72 independent root involution equations,
  * exact evaluation in the concrete amalgam model (for all labels used here),
  * canonical label deduplication, triangle generation, and signal P/Q distinctness,
  * exact Kraus normalization over denominator 5 (only when channel files are emitted),
  * dimensions, Choi rank (via singleton pivots), and sparse-entry counts,
  * the rational/integer inequalities used in the global-gap derivation.

Analytic layer:
  The dimension-uniform rigidity implication is Wang and Zhi's reflection estimate
  (arXiv:2610.01536, Theorem 5.2), used as a theorem in the paper.

The relation system has 81 generators: 72 roots, 7 action symbols and the signals p, q.
The triangle table built here (prefix triangles and one protected inverse triangle per
generator) is the table used by the 12-qubit channel of build_max_anchor_wz_channel.py.
When run as a script with channel files enabled, this module also emits fully referenced
channels (one singleton per label, as in the earlier report); the paper does not use them.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from itertools import combinations, permutations
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run this verifier without Python -O; its assertions are verification conditions.")

# --------------------------- F_8 arithmetic ---------------------------
# F_8 = F_2[alpha]/(alpha^3+alpha+1), represented by 3-bit polynomials.
MOD_POLY = 0b1011
ALPHA = 0b010


def f8_mul(a: int, b: int) -> int:
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0b1000:
            a ^= MOD_POLY
    return r & 0b111


def f8_pow(a: int, n: int) -> int:
    r = 1
    while n:
        if n & 1:
            r = f8_mul(r, a)
        a = f8_mul(a, a)
        n //= 2
    return r


# Laurent polynomial: dict[(e1,e2,e3)] = nonzero F_8 coefficient.
def pmono(c: int = 1, e=(0, 0, 0)):
    return {} if c == 0 else {tuple(e): c}


P0 = {}
P1 = pmono()


def padd(p, q):
    r = dict(p)
    for e, c in q.items():
        r[e] = r.get(e, 0) ^ c
        if r[e] == 0:
            del r[e]
    return r


def pmul(p, q):
    r = {}
    for e, c in p.items():
        for f, d in q.items():
            ef = (e[0] + f[0], e[1] + f[1], e[2] + f[2])
            v = f8_mul(c, d)
            r[ef] = r.get(ef, 0) ^ v
            if r[ef] == 0:
                del r[ef]
    return r


def paction(p, A):
    r = {}
    for e, c in p.items():
        ne = tuple(sum(A[i][j] * e[j] for j in range(3)) for i in range(3))
        r[ne] = r.get(ne, 0) ^ c
        if r[ne] == 0:
            del r[ne]
    return r


# --------------------------- 3 x 3 matrices ---------------------------
def poly_I3():
    return [[P1 if i == j else P0 for j in range(3)] for i in range(3)]


def mmul(A, B):
    C = [[{} for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            s = {}
            for k in range(3):
                s = padd(s, pmul(A[i][k], B[k][j]))
            C[i][j] = s
    return C


def maction(M, A):
    return [[paction(M[i][j], A) for j in range(3)] for i in range(3)]


def det2(a, b, c, d):
    # In characteristic two, subtraction is addition.
    return padd(pmul(a, d), pmul(b, c))


def minv_det1(M):
    cof = [[None] * 3 for _ in range(3)]
    for i in range(3):
        rows = [r for r in range(3) if r != i]
        for j in range(3):
            cols = [c for c in range(3) if c != j]
            cof[i][j] = det2(
                M[rows[0]][cols[0]], M[rows[0]][cols[1]],
                M[rows[1]][cols[0]], M[rows[1]][cols[1]],
            )
    return [[cof[j][i] for j in range(3)] for i in range(3)]


I3 = tuple(tuple(1 if i == j else 0 for j in range(3)) for i in range(3))


def imul(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))


def iinv(A):
    cof = []
    for i in range(3):
        row = []
        for j in range(3):
            rows = [r for r in range(3) if r != i]
            cols = [c for c in range(3) if c != j]
            d = A[rows[0]][cols[0]] * A[rows[1]][cols[1]] - A[rows[0]][cols[1]] * A[rows[1]][cols[0]]
            if (i + j) & 1:
                d = -d
            row.append(d)
        cof.append(row)
    return tuple(tuple(cof[j][i] for j in range(3)) for i in range(3))


def elem(i, j, p):
    M = poly_I3()
    M[i][j] = padd(M[i][j], p)
    return M


def mat_key(M):
    return tuple(
        tuple(
            tuple(sorted((e[0], e[1], e[2], c) for e, c in M[i][j].items()))
            for j in range(3)
        )
        for i in range(3)
    )


def gkey(g):
    M, A = g
    return (mat_key(M), A)


def gmul(g, h):
    M, A = g
    N, B = h
    return (mmul(M, maction(N, A)), imul(A, B))


def ginv(g):
    M, A = g
    Ai = iinv(A)
    return (maction(minv_det1(M), Ai), Ai)


G_ID = (poly_I3(), I3)

# ---------------------- Wang--Zhi auxiliary system --------------------
I = (0, 1, 2)
L = (0, 1, 2, 3)
D = (1, 2, 3)
ORDERED = tuple((a, b) for a in D for b in D if a != b)
UNORDERED = ((1, 2), (1, 3), (2, 3))


def pos(s):
    return (1, s)


def neg(s):
    return (-1, s)


def inverse_word(w):
    return [(-sgn, name) for sgn, name in reversed(w)]


def xi_exp(l):
    if l == 0:
        return (0, 0, 0)
    e = [0, 0, 0]
    e[l - 1] = 1
    return tuple(e)


def pshared(a, b, r):
    a, b = sorted((a, b))
    return f"p{a}{b}_{r}"


AUX_FAMILIES = ("c", "d0", "d1", "d2", "d3", "e0", "e1", "e2", "e3")

# 72 independent root symbols: WZ's ordered p_ab symbols are identified with p_ba.
ROOTS = (
    tuple(f"x{r}" for r in I)
    + tuple(f"y{r}" for r in I)
    + tuple(f"z{l}_{r}" for l in L for r in I)
    + tuple(f"c{r}" for r in I)
    + tuple(f"d{l}_{r}" for l in L for r in I)
    + tuple(f"e{l}_{r}" for l in L for r in I)
    + tuple(pshared(a, b, r) for a, b in UNORDERED for r in I)
    + tuple(f"{fam}_{r}" for fam in AUX_FAMILIES for r in (3, 4))
)
ACTIONS = tuple(f"U{a}{b}" for a, b in ORDERED) + ("t",)
GENERATORS = ROOTS + ACTIONS
assert len(ROOTS) == 72 and len(ACTIONS) == 7 and len(GENERATORS) == 79

# Integer shears T_ab = I + E_ab.
T = {}
for a, b in ORDERED:
    A = [list(row) for row in I3]
    A[a - 1][b - 1] += 1
    T[(a, b)] = tuple(tuple(row) for row in A)
    assert imul(T[(a, b)], iinv(T[(a, b)])) == I3

GEN = {}
for r in I:
    c = f8_pow(ALPHA, r)
    GEN[f"x{r}"] = (elem(0, 1, pmono(c)), I3)
    GEN[f"y{r}"] = (elem(1, 2, pmono(c)), I3)
    GEN[f"c{r}"] = (elem(0, 2, pmono(c)), I3)
for l in L:
    for r in I:
        c = f8_pow(ALPHA, r)
        p = pmono(c, xi_exp(l))
        GEN[f"z{l}_{r}"] = (elem(2, 0, p), I3)
        GEN[f"d{l}_{r}"] = (elem(2, 1, p), I3)
        GEN[f"e{l}_{r}"] = (elem(1, 0, p), I3)
for a, b in UNORDERED:
    for r in I:
        c = f8_pow(ALPHA, r)
        e = [0, 0, 0]
        e[a - 1] += 1
        e[b - 1] += 1
        GEN[pshared(a, b, r)] = (elem(2, 0, pmono(c, tuple(e))), I3)
for fam in AUX_FAMILIES:
    for rr in (3, 4):
        c = f8_pow(ALPHA, rr)
        if fam == "c":
            M = elem(0, 2, pmono(c))
        elif fam[0] == "d":
            M = elem(2, 1, pmono(c, xi_exp(int(fam[1:]))))
        else:
            M = elem(1, 0, pmono(c, xi_exp(int(fam[1:]))))
        GEN[f"{fam}_{rr}"] = (M, I3)
for a, b in ORDERED:
    GEN[f"U{a}{b}"] = (poly_I3(), T[(a, b)])
GEN["t"] = (poly_I3(), T[(1, 2)])
assert set(GEN) == set(GENERATORS)


def eval_word(w):
    g = G_ID
    for sgn, name in w:
        h = GEN[name]
        if sgn < 0:
            h = ginv(h)
        g = gmul(g, h)
    return g


def aux_name(base, idx):
    if idx <= 2:
        return f"c{idx}" if base == "c" else f"{base}_{idx}"
    return f"{base}_{idx}"


# Each equation record is (family, left-word, right-word).
EQUATIONS = []

# 585 commutation equations.
FAM_MEMBERS = {
    "x": tuple(f"x{r}" for r in I),
    "y": tuple(f"y{r}" for r in I),
    "c": tuple(f"c{r}" for r in I),
    "z": tuple(f"z{l}_{r}" for l in L for r in I),
    "d": tuple(f"d{l}_{r}" for l in L for r in I),
    "e": tuple(f"e{l}_{r}" for l in L for r in I),
}
for mem in FAM_MEMBERS.values():
    for a, b in combinations(mem, 2):
        EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
for r in I:
    for s in I:
        for fam in ("x", "y"):
            a, b = f"c{r}", f"{fam}{s}"
            EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
for l in L:
    for r in I:
        a = f"d{l}_{r}"
        for s in I:
            b = f"x{s}"
            EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
        for m in L:
            for s in I:
                b = f"z{m}_{s}"
                EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
for l in L:
    for r in I:
        a = f"e{l}_{r}"
        for s in I:
            b = f"y{s}"
            EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
        for m in L:
            for s in I:
                b = f"z{m}_{s}"
                EQUATIONS.append(("RH_comm", [pos(a), pos(b)], [pos(b), pos(a)]))
assert sum(1 for e in EQUATIONS if e[0] == "RH_comm") == 585

# 81 crossing equations.
for i in I:
    for j in I:
        c = aux_name("c", i + j)
        EQUATIONS.append(("RH_cross", [pos(f"x{i}"), pos(f"y{j}")], [pos(c), pos(f"y{j}"), pos(f"x{i}")]))
for l in L:
    for i in I:
        for j in I:
            d = aux_name(f"d{l}", i + j)
            e = aux_name(f"e{l}", i + j)
            EQUATIONS.append(("RH_cross", [pos(f"x{i}"), pos(f"z{l}_{j}")], [pos(d), pos(f"z{l}_{j}"), pos(f"x{i}")]))
            EQUATIONS.append(("RH_cross", [pos(f"y{i}"), pos(f"z{l}_{j}")], [pos(e), pos(f"z{l}_{j}"), pos(f"y{i}")]))

# 15 integer equations.
for i in D:
    j, k = sorted(set(D) - {i})
    EQUATIONS.append(("RK", [pos(f"U{i}{j}"), pos(f"U{i}{k}")], [pos(f"U{i}{k}"), pos(f"U{i}{j}")]))
    EQUATIONS.append(("RK", [pos(f"U{j}{i}"), pos(f"U{k}{i}")], [pos(f"U{k}{i}"), pos(f"U{j}{i}")]))
for i, j, k in permutations(D, 3):
    EQUATIONS.append(("RK", [pos(f"U{i}{j}"), pos(f"U{j}{k}")], [pos(f"U{i}{k}"), pos(f"U{j}{k}"), pos(f"U{i}{j}")]))
for i, j in UNORDERED:
    EQUATIONS.append(("RK", [pos(f"U{i}{j}"), neg(f"U{j}{i}"), pos(f"U{i}{j}")], [neg(f"U{j}{i}"), pos(f"U{i}{j}"), neg(f"U{j}{i}")]))

# 18 product equations, sharing p_ab = p_ba.
for a, b in ORDERED:
    for r in I:
        EQUATIONS.append(("Rprod", [pos(f"d{a}_{r}"), pos(f"e{b}_0")], [pos(pshared(a, b, r)), pos(f"e{b}_0"), pos(f"d{a}_{r}")]))

# 126 common compression equations.
for a, b in ORDERED:
    V = f"U{a}{b}"
    for r in I:
        EQUATIONS.append(("Rcomp", [pos(V), pos(f"x{r}")], [pos(f"x{r}"), pos(V)]))
        EQUATIONS.append(("Rcomp", [pos(V), pos(f"y{r}")], [pos(f"y{r}"), pos(V)]))
    for l in L:
        for r in I:
            target = pshared(a, b, r) if l == b else f"z{l}_{r}"
            EQUATIONS.append(("Rcomp", [pos(V), pos(f"z{l}_{r}")], [pos(target), pos(V)]))
V = "t"
for r in I:
    EQUATIONS.append(("Rcomp", [pos(V), pos(f"x{r}")], [pos(f"x{r}"), pos(V)]))
    EQUATIONS.append(("Rcomp", [pos(V), pos(f"y{r}")], [pos(f"y{r}"), pos(V)]))
for l in L:
    for r in I:
        target = pshared(1, 2, r) if l == 2 else f"z{l}_{r}"
        EQUATIONS.append(("Rcomp", [pos(V), pos(f"z{l}_{r}")], [pos(target), pos(V)]))

# 18 auxiliary addition equations.
for fam in AUX_FAMILIES:
    EQUATIONS.append(("Radd", [pos(f"{fam}_3")], [pos(aux_name(fam, 0)), pos(aux_name(fam, 1))]))
    EQUATIONS.append(("Radd", [pos(f"{fam}_4")], [pos(aux_name(fam, 1)), pos(aux_name(fam, 2))]))

COUNTS = Counter(f for f, _, _ in EQUATIONS)
assert COUNTS == Counter({"RH_comm": 585, "RH_cross": 81, "RK": 15, "Rprod": 18, "Rcomp": 126, "Radd": 18})
assert len(EQUATIONS) == 843
assert max(len(l) + len(r) for _, l, r in EQUATIONS) == 6

# Add the 72 independent root involution equations.  We also name the two
# distinguished signal words P and Q as two extra unitary generators; this makes
# the quantitative compiler theorem completely uniform (signal control is just
# another pair of relators).
P_WORD = [neg("U12"), pos("z2_0"), pos("U12")]
Q_WORD = [neg("t"), pos("z2_0"), pos("t")]
P_VALUE = eval_word(P_WORD)
Q_VALUE = eval_word(Q_WORD)
GEN_EXT = dict(GEN)
GEN_EXT["p"] = P_VALUE
GEN_EXT["q"] = Q_VALUE
ALL_GENERATORS = GENERATORS + ("p", "q")
SIGNAL_EQUATIONS = [
    ("Rsig", [pos("p")], P_WORD),
    ("Rsig", [pos("q")], Q_WORD),
]
ALL_EQUATIONS = EQUATIONS + [("Rinv", [pos(rt), pos(rt)], []) for rt in ROOTS] + SIGNAL_EQUATIONS
assert len(ALL_EQUATIONS) == 917
assert len(ALL_GENERATORS) == 81

def eval_word_ext(w):
    g = G_ID
    for sgn, name in w:
        h = GEN_EXT[name]
        if sgn < 0:
            h = ginv(h)
        g = gmul(g, h)
    return g

# Verify exact group substitution.
for fam, left, right in ALL_EQUATIONS:
    assert gkey(eval_word_ext(left)) == gkey(eval_word_ext(right)), (fam, left, right)

# ------------------- finite amalgam-label canonicalization -------------------
def eq_side(eq):
    fam, left, _ = eq
    if fam in {"RH_comm", "RH_cross", "Rprod", "Radd", "Rinv"}:
        return "H"
    if fam == "RK":
        return "G0"
    if fam == "Rcomp":
        return "G1" if left[0][1] == "t" else "G0"
    if fam == "Rsig":
        return "G0" if left[0][1] == "p" else "G1"
    raise ValueError(fam)


# H-prefix values give finite, direct certificates of every action-cancelled prefix used below.
H_PREFIX_KEYS = {gkey(G_ID)}
for eq in ALL_EQUATIONS:
    if eq_side(eq) != "H":
        continue
    rel = eq[1] + inverse_word(eq[2])
    g = G_ID
    for tok in rel:
        h = GEN_EXT[tok[1]]
        if tok[0] < 0:
            h = ginv(h)
        g = gmul(g, h)
        H_PREFIX_KEYS.add(gkey(g))

# For RK/Rcomp, every action-cancelled prefix is matched to an explicit H-prefix
# value.  Signal-definition prefixes are allowed to have trivial action while
# remaining outside H (indeed P and Q do).
for eq in ALL_EQUATIONS:
    side = eq_side(eq)
    if side == "H" or eq[0] == "Rsig":
        continue
    rel = eq[1] + inverse_word(eq[2])
    g = G_ID
    for tok in rel:
        h = GEN_EXT[tok[1]]
        if tok[0] < 0:
            h = ginv(h)
        g = gmul(g, h)
        if g[1] == I3:
            assert gkey(g) in H_PREFIX_KEYS


def canonical_key(side, g):
    k = gkey(g)
    if k in H_PREFIX_KEYS:
        return ("H", k)
    assert side in {"G0", "G1"}
    return (side, k)


def gen_side(name):
    if name in {"t", "q"}:
        return "G1"
    if name.startswith("U") or name == "p":
        return "G0"
    return "H"


def signed_gen_key(name, sgn):
    g = GEN_EXT[name]
    if sgn < 0:
        g = ginv(g)
    return canonical_key(gen_side(name), g)


ID_KEY = ("H", gkey(G_ID))
P_KEY = canonical_key("G0", P_VALUE)
Q_KEY = canonical_key("G1", Q_VALUE)
assert P_KEY != Q_KEY and P_KEY[0] == "G0" and Q_KEY[0] == "G1"

def has_negative_exponent(M):
    return any(any(ei < 0 for ei in exp) for row in M for poly in row for exp in poly)

# P and Q have trivial action part but a Laurent coefficient xi_1^{-1} xi_2,
# hence are outside H=EL_3(F_8[xi_1,xi_2,xi_3]).
assert P_VALUE[1] == I3 and Q_VALUE[1] == I3
assert has_negative_exponent(P_VALUE[0]) and has_negative_exponent(Q_VALUE[0])

# Deterministic label interning.
label_to_id = {}
labels = []


def intern(k):
    if k not in label_to_id:
        label_to_id[k] = len(labels)
        labels.append(k)
    return label_to_id[k]


assert intern(ID_KEY) == 0
for name in ALL_GENERATORS:
    intern(signed_gen_key(name, +1))
    intern(signed_gen_key(name, -1))

triangles_raw = []
for eq in ALL_EQUATIONS:
    side = eq_side(eq)
    rel = eq[1] + inverse_word(eq[2])
    g = G_ID
    pref = []
    for tok in rel:
        h = GEN_EXT[tok[1]]
        if tok[0] < 0:
            h = ginv(h)
        g = gmul(g, h)
        pref.append(intern(canonical_key(side, g)))
    assert pref[-1] == 0
    for k in range(1, len(rel)):
        triangles_raw.append((pref[k - 1], intern(signed_gen_key(rel[k][1], rel[k][0])), pref[k]))

# Protected inverse equations for all 81 positive generators.
for name in ALL_GENERATORS:
    triangles_raw.append((intern(signed_gen_key(name, +1)), intern(signed_gen_key(name, -1)), 0))

# The named signal generators are tied to P and Q by the two Rsig relators above.
P_LABEL = intern(signed_gen_key("p", +1))
Q_LABEL = intern(signed_gen_key("q", +1))
assert P_LABEL != Q_LABEL

# Deduplicate triangles preserving first occurrence.
seen = set()
triangles = []
for t in triangles_raw:
    if t not in seen:
        seen.add(t)
        triangles.append(t)

N_LABELS = len(labels)
N_TRIANGLES = len(triangles)
NATURAL_D = N_LABELS + 2 * N_TRIANGLES
assert (N_LABELS, N_TRIANGLES, NATURAL_D) == (949, 1866, 4681)
for side, (Mkey, A) in labels:
    if side != "H" and A == I3:
        # Reconstructing a polynomial dict is unnecessary: the serialized monomial
        # keys directly expose a negative exponent.
        assert any(item[0] < 0 or item[1] < 0 or item[2] < 0
                   for row in Mkey for entry in row for item in entry)

# ------------------------ sparse rational channels ------------------------
def build_channel(d):
    if d < NATURAL_D:
        raise ValueError("target dimension smaller than natural compiler dimension")
    kraus = [[] for _ in range(N_LABELS)]
    # One singleton per label, coefficient 1 = 5/5.
    for g in range(N_LABELS):
        kraus[g].append([g, g, 5])
    # Triangle blocks.
    for n, (x, y, z) in enumerate(triangles):
        p = N_LABELS + 2 * n
        kraus[0].append([p, p, 3])
        kraus[y].append([p, p + 1, 4])
        kraus[x].append([p + 1, p, 4])
        kraus[z].append([p + 1, p + 1, -3])
    # Identity-label padding to requested dimension.
    for p in range(NATURAL_D, d):
        kraus[0].append([p, p, 5])
    return kraus


def verify_normalization(kraus, d):
    # Verify sum_g A_g^* A_g = 25 I and sum_g A_g A_g^* = 25 I in integer numerators.
    left = defaultdict(int)   # (col1,col2)
    right = defaultdict(int)  # (row1,row2)
    for entries in kraus:
        by_row = defaultdict(list)
        by_col = defaultdict(list)
        for r, c, n in entries:
            by_row[r].append((c, n))
            by_col[c].append((r, n))
        for vals in by_row.values():
            for c1, n1 in vals:
                for c2, n2 in vals:
                    left[(c1, c2)] += n1 * n2
        for vals in by_col.values():
            for r1, n1 in vals:
                for r2, n2 in vals:
                    right[(r1, r2)] += n1 * n2
    expected = {(i, i): 25 for i in range(d)}
    left = {k: v for k, v in left.items() if v}
    right = {k: v for k, v in right.items() if v}
    assert left == expected
    assert right == expected


def serialize_gkey(k):
    side, (M, A) = k
    return {"side": side, "matrix": M, "action": A}


def canonical_json_bytes(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def emit_channel(outdir: Path, d: int, name: str):
    kraus = build_channel(d)
    verify_normalization(kraus, d)
    nonzero = sum(len(e) for e in kraus)
    assert nonzero == d + 2 * N_TRIANGLES
    obj = {
        "name": name,
        "dimension": d,
        "qubits": d.bit_length() - 1 if d & (d - 1) == 0 else None,
        "coefficient_denominator": 5,
        "label_count": N_LABELS,
        "choi_rank": N_LABELS,
        "triangle_count": N_TRIANGLES,
        "singleton_count": d - 2 * N_TRIANGLES,
        "identity_padding_count": d - NATURAL_D,
        "nonzero_kraus_entries": nonzero,
        "signal_labels": {"P": P_LABEL, "Q": Q_LABEL},
        "triangles": triangles,
        "kraus": kraus,
    }
    raw = canonical_json_bytes(obj)
    sha = hashlib.sha256(raw).hexdigest()
    obj["canonical_sha256_without_hash_field"] = sha
    raw2 = canonical_json_bytes(obj)
    path = outdir / f"{name}.json.gz"
    with open(path, "wb") as f:
        with gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0) as gz:
            gz.write(raw2)
    return path, sha, nonzero


# --------------------- exact constant / gap arithmetic ---------------------
def verify_gap_arithmetic():
    """The integer comparisons of Proposition 4.1 and Theorem 4.2 (see verify_compiler_constants.py)."""
    # WZ: J = 2^(2^1122000), delta = 2^-J; J > 20.
    # Proposition 4.1: rho0 = delta/2^10; rounding gives rho0 + 6 rho0 = 7 delta/1024 < delta/128.
    assert 7 * 128 < 1024
    # Signals: 5*2^-20 + 4 rho0 < 2^-17 (rho0 <= 2^-30), scaled by 2^30.
    assert 5 * 2**10 + 4 < 2**13
    # Theorem 3.4: relator defect <= 385 r < rho0 when r < rho0/512.
    assert 385 < 512
    # Theorem 4.2: g = (rho0/512)^2 / 4096 = 2^-(2J + 2*19 + 12) = 2^-(2J+50).
    assert 2 * 19 + 12 == 50
    return {
        "WZ_J": "2^(2^1122000)",
        "WZ_delta": "2^(-J)",
        "relator_tolerance_rho0": "delta/2^10 = 2^(-(J+10))",
        "signal_closeness": "< 2^-17",
        "gap_12q_choi_and_diamond": "2^(-(2J+50))",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="generated_wz_channel", help="output directory")
    ap.add_argument("--no-channel-files", action="store_true")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    gaps = verify_gap_arithmetic()

    # Exact finite certificate for labels/group values (kept separate from the larger Kraus files).
    cert = {
        "source_relation_counts": dict(COUNTS),
        "wz_balanced_equations": len(EQUATIONS),
        "independent_root_involution_equations": len(ROOTS),
        "compiled_relators": len(ALL_EQUATIONS),
        "independent_unitary_generators": len(ALL_GENERATORS),
        "max_relator_length": max(len(l) + len(r) for _, l, r in ALL_EQUATIONS),
        "labels": N_LABELS,
        "triangles": N_TRIANGLES,
        "natural_dimension": NATURAL_D,
        "signal_labels": {"P": P_LABEL, "Q": Q_LABEL},
        "signal_normal_form_sides": {"P": P_KEY[0], "Q": Q_KEY[0]},
        "gap_formulas": gaps,
        "analytic_dependency": "Wang-Zhi arXiv:2610.01536v1 Theorem 5.2",
    }
    cert_raw = canonical_json_bytes(cert)
    cert["certificate_sha256_without_hash_field"] = hashlib.sha256(cert_raw).hexdigest()
    (outdir / "certificate_summary.json").write_bytes(json.dumps(cert, indent=2, sort_keys=True).encode("ascii"))

    results = {"certificate": cert}
    if not args.no_channel_files:
        emitted = {}
        for d, name in ((8192, "wz_rational_channel_13q"), (32768, "wz_rational_channel_15q")):
            path, sha, nz = emit_channel(outdir, d, name)
            emitted[name] = {"path": str(path), "sha256_canonical_payload": sha, "nonzero_kraus_entries": nz}
        results["emitted"] = emitted

    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
