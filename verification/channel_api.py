#!/usr/bin/env python3
"""Evaluate the rational channel on sparse inputs, without dense d-by-d arrays."""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import json
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parent


class RationalChannel:
    """Exact sparse action. Coefficients may be rational or complex numbers."""
    def __init__(self, path: str | Path):
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            data = json.load(handle)
        self.dimension: int = data['dimension']
        self.rank: int = data['kraus_rank']
        self.denominator: int = data['entry_denominator']
        self._columns: dict[int, dict[int, tuple[int, int]]] = defaultdict(dict)
        for label, entries in enumerate(data['kraus']):
            for row, col, numerator in entries:
                if label in self._columns[col]:
                    raise ValueError('Expected weighted partial-permutation Kraus operators.')
                self._columns[col][label] = (row, numerator)

    def image_of_matrix_unit(self, row: int, column: int) -> dict[tuple[int, int], Fraction]:
        """Return the nonzero coefficients of Phi(|row><column|)."""
        if not (0 <= row < self.dimension and 0 <= column < self.dimension):
            raise IndexError(f'Input indices must be between 0 and {self.dimension - 1}.')
        out: dict[tuple[int, int], Fraction] = {}
        for label in self._columns[row].keys() & self._columns[column].keys():
            a, n = self._columns[row][label]
            b, k = self._columns[column][label]
            out[a, b] = out.get((a, b), Fraction(0)) + Fraction(n * k, self.denominator ** 2)
        return {index: value for index, value in out.items() if value}

    def apply_sparse(self, entries: Iterable[tuple[int, int, Fraction | complex | int]]) -> dict:
        """Apply Phi to an iterable of (row, column, coefficient) entries.

        Duplicate input entries are added. Fraction inputs give exact outputs.
        Floating-point/complex inputs keep their ordinary arithmetic semantics.
        """
        out: dict = {}
        for i, j, coefficient in entries:
            for index, value in self.image_of_matrix_unit(i, j).items():
                out[index] = out.get(index, 0) + coefficient * value
        return {index: value for index, value in out.items() if value}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--variant', choices=['one_anchor', 'two_anchor', 'fifteen_qubits'],
                   default='fifteen_qubits')
    p.add_argument('--matrix-unit', nargs=2, type=int, metavar=('ROW', 'COLUMN'), default=[2, 2])
    args = p.parse_args()
    channel = RationalChannel(ROOT / 'data' / f'{args.variant}_channel.json.gz')
    result = channel.image_of_matrix_unit(*args.matrix_unit)
    print(json.dumps({'dimension': channel.dimension, 'input_matrix_unit': args.matrix_unit,
                      'output': [[i, j, str(x)] for (i, j), x in sorted(result.items())]}, indent=2))


if __name__ == '__main__':
    main()
