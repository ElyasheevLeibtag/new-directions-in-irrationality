"""Verify the delivered source and fixed-data manifest, without dependencies."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parent

def main():
    expected=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    failures=[]
    for relative,digest in expected.items():
        path=(ROOT/relative).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            failures.append(f'Missing or invalid path: {relative}')
        elif hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            failures.append(f'Hash mismatch: {relative}')
    if failures:
        print('\n'.join(failures));raise SystemExit(1)
    print(f'PASS: {len(expected)} source and fixed-data files match the release manifest.')

if __name__=='__main__':
    main()
