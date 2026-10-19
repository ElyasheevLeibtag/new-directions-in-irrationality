"""Run each selected notebook in a separate kernel using this Python environment.

Usage: python run_notebooks.py [relative/notebook.ipynb ...] [--reverse]
No global kernelspec registration, Lean installation, or input archive is used.
Executed copies and a machine-readable report go under verification/.
"""
from pathlib import Path
import argparse
import importlib.metadata
import json
import os
import sys
import tempfile
import time

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

ROOT = Path(__file__).resolve().parent

def main():
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise RuntimeError('Assertions must remain enabled; remove -O, -OO and PYTHONOPTIMIZE.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebooks', nargs='*')
    parser.add_argument('--reverse', action='store_true', help='Run in reverse order to check independence')
    parser.add_argument('--timeout', type=int, default=3600, help='Maximum seconds per code cell')
    args = parser.parse_args()
    # Restrict discovery to delivered notebooks: a user's .venv may itself
    # contain example/test notebooks, which must never be executed here.
    delivered = list(ROOT.glob('[0-9][0-9]_*.ipynb'))
    for subject in ('joint', 'zeta2', 'pure_zeta3'):
        delivered.extend((ROOT/subject).glob('[0-9][0-9]_*.ipynb'))
    paths = [ROOT/p for p in args.notebooks] if args.notebooks else sorted(delivered)
    if args.reverse:
        paths.reverse()
    out = ROOT/'verification'
    out.mkdir(exist_ok=True)
    report = {'python': sys.version.split()[0], 'environment': {}, 'notebooks': [],
              'scope': 'Independent Python computations; no Lean build or complete formal proof.'}
    for name in ('ramanujantools','sympy','python-flint','mpmath','numpy','nbformat','nbclient','ipykernel'):
        report['environment'][name] = importlib.metadata.version(name)
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    (out/'.runtime/ipython').mkdir(parents=True, exist_ok=True)
    (out/'.runtime/jupyter').mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('IPYTHONDIR', str(out/'.runtime/ipython'))
    os.environ.setdefault('JUPYTER_RUNTIME_DIR', str(out/'.runtime/jupyter'))
    with tempfile.TemporaryDirectory(prefix='cmf-kernels-') as temp:
        kernel_dir = Path(temp)/'cmf-current'
        kernel_dir.mkdir()
        (kernel_dir/'kernel.json').write_text(json.dumps({
            'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
            'display_name': 'CMF current Python', 'language': 'python',
        }))
        for path in paths:
            path = path.resolve()
            relative = path.relative_to(ROOT)
            print(f'Running {relative}', flush=True)
            start = time.monotonic()
            nb = nbformat.read(path, as_version=4)
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    cell.outputs, cell.execution_count = [], None
            manager = KernelManager(kernel_name='cmf-current', ip='127.0.0.1',
                                    kernel_spec_manager=KernelSpecManager(kernel_dirs=[temp]))
            target = out/'executed'/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            status, error = 'passed', None
            try:
                NotebookClient(nb, km=manager, timeout=args.timeout, allow_errors=False,
                               resources={'metadata': {'path': str(path.parent)}}).execute()
            except Exception as exc:
                status, error = 'failed', str(exc)
            finally:
                if manager.has_kernel:
                    manager.shutdown_kernel(now=True)
            nbformat.write(nb, target)
            item = {'notebook': str(relative), 'status': status,
                    'seconds': round(time.monotonic()-start, 2),
                    'code_cells': sum(c.cell_type == 'code' for c in nb.cells)}
            if error:
                item['error'] = error
            report['notebooks'].append(item)
            (out/'execution_report.json').write_text(json.dumps(report, indent=2)+'\n')
            print(json.dumps(item), flush=True)
            if error:
                raise RuntimeError(f'Failed notebook: {relative}; see verification/execution_report.json')

if __name__ == '__main__':
    main()
