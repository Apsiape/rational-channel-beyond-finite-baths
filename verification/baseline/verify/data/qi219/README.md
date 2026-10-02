# QI-219 exact rational channel data

The file `channel.json.gz` contains a complete numerical quantum channel,
not a finite-dimensional control standing in for a missing hard instance.
Its nonmembership claim is CANDIDATE under the proof and external source
contract of QI-218/219. The **positive distance is not numerically evaluated**.

## Numerical map

System dimension: **44086**. Kraus count and Choi rank: **11230**.
Total nonzero Kraus entries: **76942**.

After ordinary gzip decompression, `kraus[g]` is a list of `[row,column,n]`.
The corresponding matrix entry is **n/5**, and every unlisted entry is zero.
Indices are zero-based. The map is `Phi(X)=sum_g A_g X A_g*`.
All integers `n` are in `{5,3,4,-3}`. No polynomial variable or uncomputed
trace occurs in these matrix entries.

`instance.json` includes the exact Choi spectrum, as `[weight,multiplicity]`
pairs with common denominator **1102150**. It also includes the SHA-256 of the
**uncompressed** canonical channel JSON. This avoids depending on gzip library
header conventions. A fingerprint is not an algebraic or positivity check.

## Proof data in the same object

- `labels[g]` gives a representative word in signed generator numbers and an
  exact reduced-syllable representation in the matrix amalgam.
- A factor has `side` in `{0,1}`, a row-major integer `action` matrix, and a
  row-major Laurent `matrix`. Each polynomial is a list of
  `[exponent_x1,exponent_x2,exponent_x3,coefficient_mod_5]`.
- `triangles` gives all distinct retained nontrivial group-label products.
- `generator_labels` binds each signed presentation generator to its anchor.
- The literal finite presentation is `../qi218/presentation.json`.

Reduced amalgam words do not have unique syllable lists without chosen coset
transversals. Equality is checked by reducing the product of one inverse with
the other. A polynomial factor belongs to the amalgam subgroup precisely when
its action is identity and all Laurent exponents are nonnegative.

## Reproduce

From the submission-package root:

```sh
python verify/rigid_channel_builder.py
python verify/kazhdan_compressor_presentation.py
python verify/rigid_channel_instance.py
python verify/rigid_channel_builder.py --emit rebuilt
python verification/run_checks.py
```

The builder uses only the Python standard library and contains no unimplemented
external compiler call. The independent checker uses a separate implementation
of the exact Laurent and amalgam arithmetic. Both share the mathematical
specification and Python runtime; neither verifies the external centralizer
theorem or computes a numerical separation constant.

Negative controls (each must exit nonzero):

```sh
python verify/kazhdan_compressor_presentation.py --perturb
python verify/rigid_channel_instance.py --perturb
python verify/rigid_channel_instance.py --perturb-label
python verify/rigid_channel_builder.py --perturb
```
