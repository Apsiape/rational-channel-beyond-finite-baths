# Retained baseline and reference-only reduction

This directory preserves the exact arithmetic, data, and verification suites from the supplied quantitative-progress package. It is a dependency of the revised release, not the current manuscript.

The original fully referenced channel is `verify/data/qi219/channel.json.gz`; its presentation is `verify/data/qi218/presentation.json`. The two files in `data/` are the earlier **reference-only** reductions in dimensions 32,858 and 32,857. The smaller **triangle-pruned** channels of the revised paper are instead in the release-root `data/` directory.

Run either retained suite from this directory:

```sh
python verification/run_checks.py
python verification/run_quantitative_checks.py
```

The release-root command `python verify_release.py` runs both suites automatically. Their executable code and retained exact inputs have not been altered. Execution records are regenerated on a run.

The historical report and README of the supplied package are preserved in the repository history (commit 4bc5b17, directory `provenance/`). The current paper, with its appendices, is `paper.pdf` at the repository root. The records in this directory are historical: status labels such as "CANDIDATE" in `STATUS.json`, `verification/results.json` and `verify/data/qi219/` predate the paper and are superseded by its statements. Consult `../analytical_status.json` for the scope of the analytical claims: no numerical global gap or evaluated relator tolerance is supplied.
