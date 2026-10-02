#!/usr/bin/env python3
"""Check retained release files against manifest.json (standard library only).

This detects byte changes, not mathematical errors. Verification reports are
intentionally outside the manifest. Run with --write after an intended change to
regenerate the manifest from the files in this directory.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent

EXCLUDED_DIRS = {'__pycache__', 'reports'}
EXCLUDED_FILES = {'manifest.json', 'baseline/verification/results.json',
                  'baseline/verification/quantitative_results.json'}

def write_manifest() -> int:
    files = []
    for path in sorted(ROOT.rglob('*')):
        rel = path.relative_to(ROOT).as_posix()
        parts = path.relative_to(ROOT).parts
        if (not path.is_file() or rel in EXCLUDED_FILES or EXCLUDED_DIRS & set(parts)
                or parts[0].startswith('rebuilt')):
            continue
        raw = path.read_bytes()
        files.append({'path': rel, 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()})
    manifest = {'format': 'rational-channel-release-sha256-v1', 'files': files}
    (ROOT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', newline='\n')
    print(f'WROTE: manifest of {len(files)} files')
    return 0

def main() -> int:
    if '--write' in sys.argv[1:]:
        return write_manifest()
    try:
        manifest = json.loads((ROOT/'manifest.json').read_text())
        entries = manifest['files']
        if not entries:
            raise ValueError('empty manifest')
        failures = []
        for item in entries:
            path = ROOT/item['path']
            if not path.resolve().is_relative_to(ROOT):
                raise ValueError('manifest path leaves the release directory')
            try:
                raw = path.read_bytes()
                if len(raw) != item['bytes'] or hashlib.sha256(raw).hexdigest() != item['sha256']:
                    failures.append(item['path']+': bytes differ')
            except OSError:
                failures.append(item['path']+': missing or unreadable')
        if failures:
            for failure in failures:
                print('FAIL:', failure, file=sys.stderr)
            return 1
        print(f'PASS: {len(entries)} retained files match the SHA-256 manifest')
        print('SCOPE: byte integrity only; run verify_release.py for mathematical finite-data checks')
        return 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print('FAIL: '+str(exc), file=sys.stderr)
        return 1

if __name__ == '__main__':
    sys.exit(main())
