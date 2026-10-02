#!/usr/bin/env python3
"""Replay exact instance checks, corruption controls, and deterministic emission.

These tests concern finite algebraic data. They do not verify the external
rigidity theorems or certify a numerical distance from finite tracial baths.
Only Python's standard library is required.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = (
    "kazhdan_compressor_presentation",
    "rigid_channel_builder",
    "rigid_channel_instance",
)
CONTROLS = (
    ("kazhdan_compressor_presentation", "--perturb"),
    ("rigid_channel_builder", "--perturb"),
    ("rigid_channel_instance", "--perturb"),
    ("rigid_channel_instance", "--perturb-label"),
)


def execute(arguments: list[str], timeout: int) -> dict:
    start = time.monotonic()
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    try:
        result = subprocess.run(
            [sys.executable, *arguments], cwd=ROOT, env=env,
            capture_output=True, text=True, timeout=timeout, check=False,
        )
        return {
            "command": ["python", *arguments],
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "elapsed_seconds": round(time.monotonic() - start, 3),
        }
    except subprocess.TimeoutExpired:
        return {
            "command": ["python", *arguments], "returncode": None,
            "stdout": "", "stderr": f"Timed out after {timeout} seconds.",
            "elapsed_seconds": round(time.monotonic() - start, 3),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path,
                        default=ROOT / "verification" / "results.json")
    parser.add_argument("--timeout", type=int, default=180,
                        help="Seconds allowed for each child program (default: 180).")
    parser.add_argument("--skip-reemit", action="store_true",
                        help="Skip optional temporary-directory emission comparison.")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")

    cases = []
    for name in SCRIPTS:
        record = execute([f"verify/{name}.py"], args.timeout)
        expected = (ROOT / "verify" / "expected" / f"{name}.txt").read_text()
        record.update(
            name=name, kind="positive",
            passed=(record["returncode"] == 0 and record["stdout"] == expected),
        )
        cases.append(record)
        print(("PASS" if record["passed"] else "FAIL") + f" positive {name}", flush=True)

    for optimized in (False, True):
        for name, flag in CONTROLS:
            arguments = (["-O"] if optimized else []) + [f"verify/{name}.py", flag]
            record = execute(arguments, args.timeout)
            message = record["stdout"] + record["stderr"]
            record.update(
                name=f"{name} {flag}" + (" [-O]" if optimized else ""),
                kind="negative_control", optimized=optimized,
                passed=(record["returncode"] not in (0, None) and "FAIL:" in message),
            )
            cases.append(record)
            print(("PASS" if record["passed"] else "FAIL") +
                  f" expected rejection {record['name']}", flush=True)

    if not args.skip_reemit:
        with tempfile.TemporaryDirectory(prefix="rational-channel-rebuild-") as directory:
            out = Path(directory)
            record = execute(["verify/rigid_channel_builder.py", "--emit", str(out)],
                             args.timeout)
            # Do not retain local temporary-directory names in the public report.
            record["command"][-1] = "<temporary-directory>"
            comparisons = {}
            if record["returncode"] == 0:
                paths = {
                    "presentation.json": ROOT / "verify/data/qi218/presentation.json",
                    "channel.json.gz": ROOT / "verify/data/qi219/channel.json.gz",
                    "instance.json": ROOT / "verify/data/qi219/instance.json",
                }
                for name, reference in paths.items():
                    candidate, retained = (out / name).read_bytes(), reference.read_bytes()
                    comparisons[name] = {"byte_identical": candidate == retained}
                    if name.endswith(".gz"):
                        candidate, retained = gzip.decompress(candidate), gzip.decompress(retained)
                        comparisons[name]["uncompressed_identical"] = candidate == retained
                    comparisons[name]["content_sha256"] = hashlib.sha256(candidate).hexdigest()
            # Gzip streams can differ between zlib/Python releases; uncompressed
            # canonical bytes are the mathematical binding. Also report exact
            # compressed-byte identity, but do not require it across platforms.
            canonical_ok = bool(comparisons) and all(
                item.get("uncompressed_identical", item["byte_identical"])
                for item in comparisons.values()
            )
            record.update(name="deterministic_reemission", kind="reemission",
                          comparisons=comparisons,
                          passed=(record["returncode"] == 0 and canonical_ok))
            cases.append(record)
            print(("PASS" if record["passed"] else "FAIL") + " deterministic re-emission",
                  flush=True)

    passed = sum(bool(case["passed"]) for case in cases)
    report = {
        "format": "rational-channel-letter-verification-v1",
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version.split()[0],
        "counts": {"passed": passed, "total": len(cases)},
        "scope": "Exact finite presentation, quotient word table, rational channel data, "
                 "specified corruption controls, and deterministic re-emission only.",
        "analytical_claim_status": "CANDIDATE",
        "external_rigidity_theorems_formalized": False,
        "numerical_gap_lower_bound": None,
        "cases": cases,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", newline="\n")
    print(f"RESULT: {passed}/{len(cases)} checks passed.")
    print("SCOPE: finite data only; analytical claims remain CANDIDATE; numerical gap unevaluated.")
    return 0 if passed == len(cases) else 1


if __name__ == "__main__":
    raise SystemExit(main())
