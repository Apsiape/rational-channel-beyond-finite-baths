#!/usr/bin/env python3
"""Rebuild the full presentation-to-channel chain in a separate directory.

Python 3.10+, standard library only. Existing output directories are rejected
unless --overwrite is supplied. This never modifies retained release data.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'rebuilt')
    parser.add_argument('--overwrite', action='store_true',
                        help='Allow rebuilding files in an existing output directory.')
    args = parser.parse_args()
    out = args.out.expanduser().resolve()
    protected = [ROOT/'baseline', ROOT/'data', ROOT/'manuscript', ROOT/'tests',
                 ROOT/'reports', ROOT/'provenance']
    if out == ROOT or ROOT.is_relative_to(out) or any(
            out == p or out.is_relative_to(p) for p in protected):
        parser.error('choose a separate output directory, not a retained release directory')
    if out.exists() and (not out.is_dir() or not args.overwrite):
        parser.error('output already exists; choose a new directory or use --overwrite')
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    steps = [
        [str(ROOT/'baseline/verify/rigid_channel_builder.py'), '--emit', str(out/'baseline')],
        [str(ROOT/'build_channels.py'), '--source', str(out/'baseline/channel.json.gz'),
         '--out', str(out/'data')],
        [str(ROOT/'verify_reduced.py'), '--data', str(out/'data'),
         '--json', str(out/'finite_instances.json')],
    ]
    try:
        for step in steps:
            print('RUN:', Path(step[0]).name, flush=True)
            subprocess.run([sys.executable, *step], cwd=ROOT, env=env, check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        print('FAIL: reconstruction stopped: '+str(exc), file=sys.stderr)
        return 1
    print('PASS: full reconstruction and independent reduced-data checks completed')
    print('OUTPUT:', out)
    print('SCOPE: finite data; no evaluated relator tolerance or numerical global gap')
    return 0

if __name__ == '__main__':
    sys.exit(main())
