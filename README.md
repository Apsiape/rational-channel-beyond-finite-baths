# A rational quantum channel beyond finite tracial baths

Nidhal Mghirbi and Seth Douglas — October 2026.

[Read the paper](paper.pdf) · [TeX source](paper.tex) ·
[Verification scope](verification/README.md)

A quantum channel on a finite system always has a finite realization with a pure
environment. With a maximally mixed environment this can fail even approximately: since
MIP* = RE, some factorizable channels lie outside the closure of the finite matrix-tracial
factorizations (Haagerup and Musat). This paper constructs an explicit one: a unital channel
on fifteen qubits with 11,230 Kraus operators, whose entries are integers divided by five,
that factorizes exactly through a type II<sub>1</sub> factor but cannot be approximated with
any finite maximally mixed environment.

The construction starts from an explicit finite presentation, with 43 generators and 4,417
relators, of a group that maps into the amalgamated double of a Kun–Thom pair over the field
with five elements. In it, Thom's witness commutator becomes a word of length six, which every
tracial matrix-ultraproduct representation kills by the centralizer-rigidity theorems of
Alekseev, Liu and Thom and of Thom (recent preprints). Rational two-level blocks enforce the
relations, and the Choi matrix gives the word a trace the identity cannot have.

The paper also proves two quantitative statements:

- Against finite tracial channels whose Choi support lies inside the target's, the
  fifteen-qubit channel is at normalized Choi distance more than 1/316 and half-diamond
  distance at least 1/2.
- For the original channel, which keeps a reference position for every label, the distance
  Δ_J to the closure is positive, and serving n uses in the causal model of the companion
  papers needs memory at least nΔ_J²/(24 ln 2) − 2 at error at most min(1/16, Δ_J/2) and with
  at most (n/2)·log₂ 11,230 exchanged qubits. A uniform relator tolerance, which exists but is
  not evaluated, would give an explicit lower bound on Δ_J.

The appendices give the full presentation, the word and channel algorithms, all proofs and the
data format. This repository contains the exact data, programs that rebuild every channel from
the presentation, and independent finite checks of the scope stated in
[verification/README.md](verification/README.md).

## Reproduce

The finite checks need Python 3.10 or later and its standard library only. From
`verification/`:

```sh
python verify_manifest.py     # byte integrity of the retained files
python verify_release.py      # all finite-data checks, including the two baseline suites
python rebuild_all.py --out rebuilt
```

The second command runs 25 top-level checks, two of which run nested suites, including
deliberate-corruption controls in ordinary and optimized Python; it takes about a minute.
The third regenerates the original presentation and channel and every reduced channel, then
checks them independently. See [verification/README.md](verification/README.md).

`python build.py` rebuilds the paper, with its appendices, with pdflatex.

## License and citation

The manuscript and repository content are available under **CC BY 4.0**.
The software and machine-readable data are additionally available under
the **MIT License**, at your option; see [RIGHTS.md](RIGHTS.md) for the scope.

Please cite Nidhal Mghirbi and Seth Douglas, *A rational quantum channel beyond finite
tracial baths* (2026), version 1.0. [CITATION.cff](CITATION.cff) provides
machine-readable citation metadata. A DOI will be added when the first version is archived.
