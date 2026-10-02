# Verification

Everything here runs with Python 3.10 or later and its standard library only. No network,
optimizer or third-party package is used. Run the commands from this folder.

## What is checked

```sh
python verify_manifest.py
python verify_release.py
```

`verify_manifest.py` checks the retained files against `manifest.json` (byte integrity, not
mathematics). After an intended change, `python verify_manifest.py --write` regenerates the
manifest.

`verify_release.py` runs 25 top-level checks and writes its record to
`reports/verification.json`. Two of the checks run nested suites (12 baseline checks and 17
quantitative checks). The counts include deliberate-corruption controls, which pass when altered
inputs are rejected; they are counts of tests, not of theorems. The suite checks:

- the original exact presentation and the quotient arithmetic;
- the canonical label table and all 438 free-group deletion certificates, each built from
  retained equations only, and the protected inverse equations;
- label coverage and an independent reconstruction of the three reduced channels;
- both Kraus normalization identities and the complete nonzero Choi spectrum;
- the exact scalar inequalities behind the support-restricted bounds;
- reproducible emission, and rejection of eight kinds of corrupted data in ordinary and
  optimized Python.

`python verify_release.py --quick --report reports/quick_verification.json` omits the two
baseline suites and is not the complete verification.

## Rebuild every channel from the finite presentation

```sh
python rebuild_all.py --out rebuilt
python verify_reduced.py --data rebuilt/data
```

The first command regenerates the original presentation and channel, emits the reduced
channels and runs the independent reduced-data checker; it refuses an existing output folder
unless `--overwrite` is given. The emitter and the checker share no pruning or
matrix-construction code: the checker derives the reduced matrices by deleting and remapping
slots of the original sparse channel. Re-emission is compared on canonical uncompressed bytes,
since gzip streams vary across library versions.

## Use a channel without dense matrices

```sh
python channel_api.py --variant fifteen_qubits --matrix-unit 0 1
python channel_api.py --variant fifteen_qubits --matrix-unit 2 2
```

The first result is zero; the second evaluates a nontrivial triangle column. In Python:

```python
from pathlib import Path
from fractions import Fraction
from channel_api import RationalChannel

channel = RationalChannel(Path("data/fifteen_qubits_channel.json.gz"))
image = channel.image_of_matrix_unit(2, 2)   # maps (row, column) to a Fraction
```

Inputs with exact `Fraction` values are evaluated exactly.

## Files

- `data/`: the three reduced channel tables, all deletion certificates, and construction
  summaries. A sparse entry `[row, column, numerator]` in `kraus[g]` is that matrix entry divided
  by five; indices are zero-based. Group generator 43 names the witness word; its channel label
  is 85. Label 0 is the identity.
- `baseline/`: the original presentation, the original fully referenced channel, its builder
  and an independent checker, and the earlier reference-only reductions.
- `build_channels.py`, `verify_reduced.py`, `tests/exact_lemmas.py`: the reduced-channel emitter,
  the independent reduced-data checker, and exact scalar and probe checks.
- `analytical_status.json`: which statements are finite checks, which rest on external theorems,
  and which are conditional; the global numerical gap and the relator tolerance are recorded as
  not evaluated.

## Limits

The code does not formalize property (T), the external rigidity theorems, the ultraproduct
extraction or the analytic inequalities; those arguments are written out or attributed in the
paper and its appendices. Passing finite tests does not certify them. No minimal dimension, evaluated
relator tolerance or global robustness constant is claimed.
