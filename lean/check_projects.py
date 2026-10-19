"""Build the pinned Lean projects and print the exported theorems' axioms.

This driver is packaging infrastructure. Its presence is not evidence that
the projects compile. Use --plan to inspect commands without changing files.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import tomllib

ROOT = Path(__file__).resolve().parents[1]
PROJECTS = [
    ('meijerg_zeta3_zeta4', 'lean', [
        ('Meijer5255/Blueprint_5255.lean', ['Meijer5255.main_joint', 'Meijer5255.main_pure3'])]),
    ('zeta2_four_trajectory_matrices', '', [
        ('Zeta2FourTrajectoryMatrices/MainTheorem.lean', ['zeta2_irrationality_measure'])]),
    ('zetas_irr', '', [
        ('blueprint.lean', ['ZetasIrr.theorem_1_6', 'ZetasIrr.irrational_combination']),
        ('pfq54_formalization_handoff/Blueprint_pfq54.lean', ['Pfq54.main_theorem'])]),
]

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', action='store_true', help='print the build plan without executing it')
    ap.add_argument('--project', choices=[p[0] for p in PROJECTS], help='check just one project')
    ap.add_argument('--timeout', type=int, default=7200, help='timeout in seconds for each command')
    args = ap.parse_args()
    if args.timeout <= 0:
        ap.error('--timeout must be positive')
    jobs = []
    for name, subdir, entrypoints in PROJECTS:
        if args.project and name != args.project:
            continue
        project = ROOT / 'lean' / name / subdir
        config = tomllib.loads((project / 'lakefile.toml').read_text())
        targets = [lib['name'] for lib in config['lean_lib']]
        jobs.append((name, project, targets, entrypoints))
    if args.plan:
        for name, project, targets, entrypoints in jobs:
            print(f'\n{name}: {project.relative_to(ROOT)}')
            print('  lake exe cache get')
            print('  lake build ' + ' '.join(targets))
            for source, theorems in entrypoints:
                print(f'  lake env lean <temporary copy of {source} with appended #print axioms>')
                print('    axioms:', ', '.join(theorems))
        return

    out = ROOT / 'verification/lean_build'
    out.mkdir(parents=True, exist_ok=True)
    report = {'status': 'running', 'scope': 'Build and axiom output only. Explicit upstream PNT assumptions are retained.', 'commands': []}
    def save():
        (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    def run(name, label, project, command):
        log = out / f'{name}_{label}.log'
        start = time.monotonic()
        row = {'project': name, 'label': label, 'command': command, 'log': log.relative_to(ROOT).as_posix()}
        report['commands'].append(row)
        print(f'{name}: {label}; log: {log}', flush=True)
        try:
            with log.open('w') as handle:
                completed = subprocess.run(command, cwd=project, stdout=handle, stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
            row['exit_code'] = completed.returncode
            row['status'] = 'passed' if completed.returncode == 0 else 'failed'
            if completed.returncode:
                raise RuntimeError(f'{name}: {label} failed; see {log}')
        except (OSError, subprocess.TimeoutExpired) as exc:
            row['status'] = 'failed'
            row['error'] = str(exc)
            raise
        finally:
            row['seconds'] = round(time.monotonic() - start, 2)
            save()
    save()
    try:
        for name, project, targets, entrypoints in jobs:
            manifest = project / 'lake-manifest.json'
            before = manifest.read_bytes()
            run(name, 'dependencies', project, ['lake', 'exe', 'cache', 'get'])
            if manifest.read_bytes() != before:
                raise RuntimeError(f'{name}: dependency manifest changed; inspect it before continuing')
            run(name, 'libraries', project, ['lake', 'build', *targets])
            for i, (source, theorems) in enumerate(entrypoints):
                text = (project / source).read_text() + '\n' + '\n'.join(f'#print axioms {t}' for t in theorems) + '\n'
                with tempfile.NamedTemporaryFile(mode='w', suffix='.lean', prefix='SupplementCheck_', dir=project, delete=False) as tmp:
                    tmp.write(text)
                    path = Path(tmp.name)
                try:
                    run(name, f'theorem_{i+1}', project, ['lake', 'env', 'lean', path.name])
                finally:
                    path.unlink(missing_ok=True)
        report['status'] = 'passed'
    except (Exception, KeyboardInterrupt) as exc:
        report['status'] = 'failed'
        report['error'] = str(exc) or 'Interrupted'
        print(report['error'], file=sys.stderr)
    finally:
        save()
    if report['status'] != 'passed':
        raise SystemExit(1)
    print('Build commands passed. Review the axiom logs; upstream PNT assumptions remain.')

if __name__ == '__main__':
    main()
