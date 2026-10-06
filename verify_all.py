#!/usr/bin/env python3
"""Run every finite verifier of the release and record the results.

Steps (Section 8 of the paper):
  1. relation_compiler.py rebuilds the relation system and its triangle table in exact arithmetic;
  2. build_max_anchor_wz_channel.py rebuilds the 12-qubit payload in a scratch directory, and its
     decompressed content must equal the shipped payload's (gzip output differs between zlib builds);
  3. the remaining scripts check the payload, the constants, the earlier channels' structure and the
     spectral certificates.
"""
import gzip
import hashlib
from datafile import read_data
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent
PAYLOAD = 'wz_max_anchor_channel_12q.json.gz'
records = []
TMP = None


def run(args, label):
    p = subprocess.run([sys.executable] + args, cwd=ROOT, text=True, capture_output=True)
    out, err = p.stdout, p.stderr
    if TMP:
        out, err = out.replace(TMP, '<scratch>'), err.replace(TMP, '<scratch>')
    records.append({'step': label, 'returncode': p.returncode, 'stdout': out, 'stderr': err})
    sys.stdout.write(f'===== {label} =====\n{p.stdout}')
    if p.stderr:
        sys.stderr.write(p.stderr)
    if p.returncode:
        raise SystemExit(p.returncode)


with tempfile.TemporaryDirectory() as tmp:
    TMP = tmp
    run(['relation_compiler.py', '--no-channel-files', '--out', tmp], 'relation_compiler.py')
    rebuilt = pathlib.Path(tmp) / PAYLOAD
    run(['build_max_anchor_wz_channel.py', str(rebuilt)], 'build_max_anchor_wz_channel.py')
    same = gzip.decompress(rebuilt.read_bytes()) == read_data(ROOT / PAYLOAD)
    records.append({'step': 'rebuilt payload equals shipped payload', 'returncode': 0 if same else 1})
    print('===== rebuilt payload equals shipped payload =====')
    print('OK identical decompressed content' if same else 'MISMATCH')
    if not same:
        raise SystemExit(1)

for s in ['verify_payload_table.py',
          'verify_channel_payload.py',
          'verify_compiler_constants.py',
          'verify_two_test_certificate.py',
          'verify_original15_leakage_constants.py',
          'verify_earlier_channels.py',
          'verify_nt_certificate.py',
          'verify_rigidity_constants.py',
          'verify_f5_old15_rigidity.py',
          'earlier/model/check_model.py',
          'earlier/model/check_labels.py']:
    run([s], s)

files = ['nt_sos_certificate.json.gz', 'nt_word_reductions.json.gz', PAYLOAD, 'f5_root_sos.json',
         'st3_opposite_reductions.json', 'earlier/presentation.json', 'earlier/original_channel.json.gz',
         'earlier/fifteen_qubits_channel.json.gz', 'earlier/two_anchor_channel.json.gz',
         'earlier/one_anchor_channel.json.gz', 'earlier/deletion_certificates.json']
# SHA-256 of the decompressed content of each .gz data file, used when a file arrives decompressed.
CONTENT_SHA256 = {
    'nt_sos_certificate.json.gz': '20b4a5cc1f4ac603cd16956bec0e84c944f47b9d361a4c39337ac730eb4de111',
    'nt_word_reductions.json.gz': 'beaa1fccdaef93db089e65e3a58491480284b2f44b7dd4795db7923576f4c52c',
    PAYLOAD: '454afe1fb21166e26af5a847eeb197e643a60b76353d16195bcde8347bca30ca',
    'earlier/original_channel.json.gz': '38150decf3d203e3eabbc5dee8d8a3f2d6390aa960c50f90fe4c87b804729b68',
    'earlier/fifteen_qubits_channel.json.gz': '1236ae6e7d587c738baac7c0e5815c90ece46ff8c9392fd0012614aae6e1aa1b',
    'earlier/two_anchor_channel.json.gz': 'c11f651e46b8dee8403b9b44a88b74a2670a0c1593e8d38756c22605b91cfbbc',
    'earlier/one_anchor_channel.json.gz': '3e56e422def2967be9561fdc4d4b97bf0b9bf534cdc0640bc0a525bf6c0c00d4',
}
hashes = {f: hashlib.sha256((ROOT / f).read_bytes()).hexdigest() for f in files if (ROOT / f).exists()}
decompressed = [f for f in files if not (ROOT / f).exists()]
for f in decompressed:
    if f not in CONTENT_SHA256 or hashlib.sha256(read_data(ROOT / f)).hexdigest() != CONTENT_SHA256[f]:
        print('MISMATCH: ' + f)
        raise SystemExit(1)
listed = {}
for line in (ROOT / 'SHA256SUMS.txt').read_text().splitlines():
    h, _, name = line.partition(' ')
    listed[name.lstrip('*').strip()] = h
if {k: v for k, v in listed.items() if k not in decompressed} != hashes:
    print('MISMATCH with SHA256SUMS.txt')
    raise SystemExit(1)
records.append({'step': 'data files match SHA256SUMS.txt'
                + (' (decompressed files checked by content hash)' if decompressed else ''), 'returncode': 0})
report = {'ok': True, 'python': sys.version, 'steps': records, 'sha256': hashes}
(ROOT / 'verification_report.json').write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
print('===== SHA256 =====')
for k, v in hashes.items():
    print(v, k)
print('OK ALL')
