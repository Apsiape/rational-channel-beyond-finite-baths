#!/usr/bin/env python3
"""Independently check pruned data by remapping the original sparse channel.

No builder module is imported. Deleted equations are proved by free reduction.
This checks finite data, not the external rigidity theorems or a global gap.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
BASE_HASH = '38150decf3d203e3eabbc5dee8d8a3f2d6390aa960c50f90fe4c87b804729b68'
VARIANTS = {
    'one_anchor': (31981, 63961, 0, Fraction(1, 527)),
    'two_anchor': (31982, 63962, 0, Fraction(1, 329)),
    'fifteen_qubits': (32768, 64748, 786, Fraction(1, 316)),
}


def require(value: bool, message: str) -> None:
    if not value:
        raise ValueError(message)


def reduce_word(letters: list[int]) -> list[int]:
    stack: list[int] = []
    for letter in letters:
        require(type(letter) is int and letter != 0, 'invalid free-group letter')
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return stack


def inverse(letters: list[int]) -> list[int]:
    return [-k for k in letters[::-1]]


def label_word(label: int) -> list[int]:
    return [label + 1] if label else []


def relation(triangle: list[int]) -> list[int]:
    x, y, z = triangle
    return reduce_word(label_word(x) + label_word(y) + inverse(label_word(z)))


def load_compressed(path: Path) -> tuple[dict, bytes]:
    raw = gzip.decompress(path.read_bytes())
    return json.loads(raw), raw


def check(data_dir: Path, mutation: str | None = None) -> dict:
    original, original_raw = load_compressed(ROOT / 'baseline/verify/data/qi219/channel.json.gz')
    require(hashlib.sha256(original_raw).hexdigest() == BASE_HASH, 'baseline fingerprint changed')
    triangles = original['triangles']
    labels = original['labels']
    generator_labels = dict(original['generator_labels'])
    require((len(triangles), len(labels), generator_labels[43]) == (16428, 11230, 85),
            'baseline dimensions or witness changed')
    certificates = json.loads((data_dir / 'deletion_certificates.json').read_text())
    summary = json.loads((data_dir / 'summary.json').read_text())
    keep = certificates['retained_original_indices']
    certs = certificates['certificates']
    if mutation == 'certificate':
        certs[0]['cells'][0]['sign'] *= -1
    if mutation == 'dependency':
        certs[0]['cells'][0]['reference'] = certs[0]['deleted']
    if mutation == 'scope':
        summary['variants']['two_anchor']['global_numerical_gap'] = '1/2'
    require(certificates['baseline_sha256'] == BASE_HASH, 'certificate source changed')
    require(certificates['original_triangle_count'] == len(triangles), 'certificate source size')
    require(keep == sorted(set(keep)) and len(keep) == 15990, 'retained index list')
    require(all(type(i) is int and 0 <= i < len(triangles) for i in keep), 'retained index range')
    deleted = [c['deleted'] for c in certs]
    require(len(deleted) == len(set(deleted)) == 438, 'deletion count')
    require(set(deleted) == set(range(len(triangles))) - set(keep), 'deletion complement')
    require(Counter(c['kind'] for c in certs) ==
            {'right_inverse': 408, 'reversed_inverse': 30}, 'deletion family counts')
    kept = set(keep)
    for certificate in certs:
        product: list[int] = []
        for cell in certificate['cells']:
            ref, sign, conjugator = cell['reference'], cell['sign'], cell['conjugator']
            require(ref in kept, 'certificate uses a deleted equation')
            require(sign in (-1, 1), 'invalid certificate sign')
            equation = relation(triangles[ref])
            if sign == -1:
                equation = inverse(equation)
            product += conjugator + equation + inverse(conjugator)
        require(reduce_word(product) == relation(triangles[certificate['deleted']]),
                'free-group certificate failed')
    retained = [triangles[i] for i in keep]
    retained_set = set(map(tuple, retained))
    for s in range(1, 44):
        require((generator_labels[s], generator_labels[-s], 0) in retained_set,
                'positive-generator inverse equation was deleted')
    require({g for t in retained for g in t} == set(range(len(labels))), 'unitary coverage')
    require(summary['total_deleted'] == 438 and summary['retained_triangles'] == 15990 and
            summary['certificate_count'] == 438 and summary['all_labels_covered'] is True and
            summary['right_inverse_pair_deletions'] == 408 and
            summary['reversed_generator_inverse_deletions'] == 30, 'pruning summary mismatch')
    results = {}
    for name, (dimension, nnz, padding, rational_bound) in VARIANTS.items():
        channel, raw = load_compressed(data_dir / f'{name}_channel.json.gz')
        if name == 'two_anchor':
            if mutation == 'coefficient':
                channel['kraus'][0][0][2] = 4
            if mutation == 'triangle':
                channel['triangles'].pop()
            if mutation == 'witness':
                channel['anchor_labels'][1] = 0
            if mutation == 'denominator':
                channel['entry_denominator'] = 7
        if name == 'fifteen_qubits' and mutation == 'padding':
            channel['identity_padding_count'] -= 1
        anchors = [0] if name == 'one_anchor' else [0, generator_labels[43]]
        require(channel['anchor_labels'] == anchors, 'reference labels mismatch')
        require(channel['entry_denominator'] == 5, 'coefficient denominator mismatch')
        require(channel['identity_padding_count'] == padding, 'identity padding mismatch')
        require(channel['dimension'] == dimension, 'system dimension mismatch')
        require(channel['baseline_uncompressed_sha256'] == BASE_HASH, 'channel source binding')
        require(channel['triangles'] == retained, 'retained triangle table mismatch')
        require(channel['labels'] == labels and channel['generator_labels'] == original['generator_labels'],
                'original label table changed')
        require(channel['kraus_rank'] == len(labels) and len(channel['kraus']) == len(labels),
                'Kraus count mismatch')
        # Remap original system slots. This is deliberately different from the builder.
        index_map = {label: i for i, label in enumerate(anchors)}
        for new_index, old_index in enumerate(keep):
            for delta in (0, 1):
                index_map[len(labels) + 2 * old_index + delta] = len(anchors) + 2 * new_index + delta
        expected = []
        for entries in original['kraus']:
            expected.append([[index_map[a], index_map[b], n] for a, b, n in entries
                             if a in index_map and b in index_map])
        natural_dimension = len(anchors) + 2 * len(keep)
        expected[0] += [[i, i, 5] for i in range(natural_dimension, dimension)]
        require(channel['kraus'] == expected, 'remapped Kraus matrices mismatch')
        rows = [0] * dimension
        columns = [0] * dimension
        slots: set[tuple[int, int]] = set()
        weights = []
        for entries in channel['kraus']:
            require(bool(entries), 'empty Kraus support')
            kraus_rows: set[int] = set()
            kraus_columns: set[int] = set()
            for a, b, n in entries:
                require(0 <= a < dimension and 0 <= b < dimension and n in (5, 3, 4, -3),
                        'invalid sparse entry')
                require((a, b) not in slots, 'overlapping Kraus supports')
                require(a not in kraus_rows and b not in kraus_columns, 'not a partial permutation')
                slots.add((a, b)); kraus_rows.add(a); kraus_columns.add(b)
                rows[a] += n * n; columns[b] += n * n
            weights.append(sum(n * n for _, _, n in entries))
        require(rows == [25] * dimension and columns == [25] * dimension, 'normalization identity')
        require(len(slots) == nnz and sum(weights) == 25 * dimension, 'sparse or spectral count')
        info = summary['variants'][name]
        require(info['global_numerical_gap'] is None, 'global numerical gap must remain null')
        require(info['dimension'] == dimension and info['rank'] == len(labels) and
                info['sparse_entries'] == nnz and info['identity_padding_count'] == padding,
                'variant summary mismatch')
        squared_gap = Fraction(weights[0] * weights[85], (25 * dimension) ** 2)
        require(squared_gap == Fraction(info['restricted_choi_gap_squared']), 'restricted Choi value')
        require(rational_bound ** 2 < squared_gap, 'rational restricted separator')
        require((weights[0], weights[85]) ==
                (info['identity_weight_numerator'], info['j_weight_numerator']), 'label weights')
        require(hashlib.sha256(raw).hexdigest() == info['sha256_uncompressed'], 'new channel fingerprint')
        if len(anchors) == 2:
            # These are the only entries in the two reference columns.
            for col, label in enumerate(anchors):
                entries = [(g, a, n) for g, k in enumerate(channel['kraus'])
                           for a, b, n in k if b == col]
                require(entries == [(label, col, 5)], 'two-position dephasing probe')
        hist = Counter(weights)
        results[name] = {
            'dimension': dimension, 'rank': len(labels), 'nonzero_entries': nnz,
            'choi_eigenvalue_denominator': 25 * dimension,
            'choi_weight_histogram': [[a, hist[a]] for a in sorted(hist)],
            'identity_weight_numerator': weights[0], 'witness_weight_numerator': weights[85],
            'restricted_choi_gap_lower_bound': str(rational_bound),
            'restricted_choi_gap_squared': str(squared_gap),
            'restricted_half_diamond_lower_bound': '1/2' if len(anchors) == 2 else None,
            'global_numerical_gap': None,
        }
    return {'deleted_equations': 438, 'verified_free_group_certificates': 438,
            'retained_equations': 15990, 'variants': results,
            'scope': 'Finite arithmetic and free-group implications, not analytical rigidity certification.'}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data')
    parser.add_argument('--json', type=Path, help='Write a machine-readable check summary.')
    parser.add_argument('--perturb', choices=['coefficient', 'triangle', 'certificate', 'dependency',
                                             'scope', 'witness', 'padding', 'denominator'])
    args = parser.parse_args()
    try:
        results = check(args.data, args.perturb)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(json.dumps(results, indent=2) + '\n', newline='\n')
        print('PASS: 438 free-group deletion certificates; 15990 retained equations; all labels covered')
        for name, v in results['variants'].items():
            print(f"PASS: {name}: d={v['dimension']}, rank={v['rank']}, entries={v['nonzero_entries']}, "
                  f"restricted Choi bound > {v['restricted_choi_gap_lower_bound']}")
        print('SCOPE: exact finite data; no numerical global gap; analytical proofs are not formalized')
        return 0
    except (ValueError, KeyError, TypeError, IndexError, OSError) as exc:
        print('FAIL: ' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
