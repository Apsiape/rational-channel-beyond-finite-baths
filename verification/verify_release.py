#!/usr/bin/env python3
"""Run the complete finite-data verification for this release.

Requires Python 3.10+ and its standard library. No network access is used.
By default both baseline suites, independent pruning checks, deliberate
corruption controls, and reproducible rebuilding are run. --quick omits only
the two baseline suites; it is not the complete verification.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--report', type=Path, default=ROOT/'reports/verification.json')
    parser.add_argument('--timeout', type=int, default=600,
                        help='Per-child timeout in seconds; baseline suites contain several checks.')
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error('--timeout must be positive')
    cases: list[dict] = []
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')

    def execute(name: str, arguments: list[str], *, cwd: Path = ROOT,
                rejection: str | None = None) -> dict:
        start = time.monotonic()
        try:
            p = subprocess.run([sys.executable, *arguments], cwd=cwd, env=env,
                               text=True, capture_output=True, timeout=args.timeout, check=False)
            passed = (p.returncode == 0) if rejection is None else (
                p.returncode != 0 and 'FAIL:' in p.stderr and rejection in p.stderr)
            record = {'name': name, 'arguments': arguments, 'working_directory': str(cwd.relative_to(ROOT)),
                      'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr,
                      'passed': passed, 'expected_rejection': rejection}
        except subprocess.TimeoutExpired:
            record = {'name': name, 'arguments': arguments, 'passed': False, 'reason': 'timeout'}
        record['seconds'] = round(time.monotonic()-start, 3)
        cases.append(record)
        print(('PASS ' if record['passed'] else 'FAIL ')+name, flush=True)
        return record

    if not args.quick:
        base = execute('baseline exact-data suite', ['verification/run_checks.py'], cwd=ROOT/'baseline')
        quant = execute('supplied quantitative-progress suite',
                        ['verification/run_quantitative_checks.py'], cwd=ROOT/'baseline')
        for rec, filename in ((base, 'results.json'), (quant, 'quantitative_results.json')):
            path = ROOT/'baseline/verification'/filename
            if path.exists():
                nested = json.loads(path.read_text())
                rec['nested_counts'] = nested['counts']
                rec['nested_report'] = str(path.relative_to(ROOT))
                rec['passed'] = rec['passed'] and nested['counts']['passed'] == nested['counts']['total']
    for optimized in (False, True):
        prefix = ['-O'] if optimized else []
        suffix = ' [optimized]' if optimized else ''
        execute('independent reduced-data checks'+suffix,
                prefix+['verify_reduced.py']+([] if optimized else ['--json','reports/finite_instances.json']))
        execute('scalar and probe checks'+suffix, prefix+['tests/exact_lemmas.py'])
        expected = {
            'coefficient': 'remapped Kraus', 'triangle': 'retained triangle table',
            'certificate': 'free-group certificate', 'dependency': 'deleted equation',
            'scope': 'global numerical gap', 'witness': 'reference labels',
            'padding': 'identity padding', 'denominator': 'coefficient denominator',
        }
        for mutation, message in expected.items():
            execute('reject '+mutation+suffix, prefix+['verify_reduced.py','--perturb',mutation],
                    rejection=message)
    with tempfile.TemporaryDirectory(prefix='rational-channel-release-') as tmp:
        directories = []
        for optimized in (False, True):
            out = Path(tmp)/('optimized' if optimized else 'normal')
            rec = execute('rebuild reduced channels'+(' [optimized]' if optimized else ''),
                          (['-O'] if optimized else [])+['build_channels.py','--out',str(out)])
            rec['arguments'][-1] = '<temporary-directory>'
            directories.append(out)
        bindings = {}
        for filename in ('deletion_certificates.json','summary.json','one_anchor_channel.json.gz',
                         'two_anchor_channel.json.gz','fifteen_qubits_channel.json.gz'):
            try:
                parts = [(directory/filename).read_bytes() for directory in directories+[ROOT/'data']]
                if filename.endswith('.gz'):
                    parts = [gzip.decompress(part) for part in parts]
                bindings[filename] = parts[0] == parts[1] == parts[2]
            except OSError:
                bindings[filename] = False
        cases.append({'name':'canonical re-emission comparison', 'passed':all(bindings.values()),
                      'bindings':bindings})
        print(('PASS ' if all(bindings.values()) else 'FAIL ')+'canonical re-emission comparison', flush=True)
    total = len(cases)
    passed = sum(bool(c['passed']) for c in cases)
    nested = [c['nested_counts'] for c in cases if 'nested_counts' in c]
    report = {
        'format':'rational-channel-revised-verification-v1',
        'executed_at_utc':datetime.now(timezone.utc).isoformat(),
        'python_version':sys.version.split()[0],
        'mode':'quick' if args.quick else 'complete',
        'counts':{'passed':passed,'total':total},
        'nested_baseline_counts':nested,
        'scope':'Finite algebra, exact rational channel data, free-group certificates, specified negative controls, and reproducible emission.',
        'analytical_proofs_formalized':False, 'external_rigidity_formalized':False,
        'global_numerical_gap':None, 'relator_tolerance_e0':None,
        'cases':cases,
    }
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n', newline='\n')
    print(f'RESULT: {passed}/{total} release checks passed.',flush=True)
    print('SCOPE: finite verification only; no evaluated global robustness constant.',flush=True)
    return 0 if passed == total else 1


if __name__ == '__main__':
    sys.exit(main())
