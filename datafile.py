"""Read a shipped data file whether it is stored compressed (name.json.gz) or decompressed (name.json).

Standard library only.  Some download and upload paths decompress .gz files; the verifiers accept
either form and always work with the decompressed content.
"""
import gzip
import os


def read_data(path):
    """Return the decompressed bytes of the data file `path`, whose name may end in .gz."""
    path = os.fspath(path)
    if os.path.exists(path):
        with open(path, 'rb') as f:
            raw = f.read()
        return gzip.decompress(raw) if path.endswith('.gz') else raw
    if path.endswith('.gz') and os.path.exists(path[:-3]):
        with open(path[:-3], 'rb') as f:
            return f.read()
    raise FileNotFoundError(path)
