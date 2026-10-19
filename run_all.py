"""Run the computational supplement with this Python interpreter.

Default: nine independent Section 6 notebooks and additional script checks.
Use --section6-only or --additional-only to select a component.
"""
from pathlib import Path
import argparse
import datetime
import importlib.metadata
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent
SCRIPTS=[
    'additional/analytic/verify_4f3_meijerg_coboundary.py',
    'additional/analytic/verify_5f4_meijerg_recessive.py',
    'additional/pi_integrals/pi_I_axis_certificate.py',
    'additional/pi_integrals/pi_I_initial_check.py',
    'additional/pi_integrals/pi_I_uniqueness_check.py',
    'additional/pi_integrals/pi_K_seed_check.py',
    'additional/recoveries/recovery_zeta5_readout_check.py',
    'additional/gregory/verify_gregory.py',
]

def main():
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise RuntimeError('Assertions must remain enabled: remove -O, -OO and PYTHONOPTIMIZE.')
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--section6-only',action='store_true')
    group.add_argument('--additional-only',action='store_true')
    parser.add_argument('--timeout',type=int,default=1800,help='Timeout per script or notebook-suite process, seconds')
    args=parser.parse_args()
    jobs=[]
    if not args.additional_only:
        jobs.append(('section6_notebooks',[sys.executable,str(ROOT/'python/run_notebooks.py'),'--reverse']))
    if not args.section6_only:
        jobs.extend((Path(p).stem,[sys.executable,str(ROOT/p)]) for p in SCRIPTS)
    out=ROOT/'verification';out.mkdir(exist_ok=True)
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONUNBUFFERED='1')
    report={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'python':sys.version.split()[0],'environment':{},'jobs':[],
            'expected_jobs':len(jobs),'status':'running',
            'scope':'Python checks only. Lean projects are packaged separately and are not executed by this runner.'}
    for name in ('ramanujantools','sympy','python-flint','mpmath','numpy','nbformat','nbclient','ipykernel'):
        report['environment'][name]=importlib.metadata.version(name)
    report_name='additional_run.json' if args.additional_only else 'section6_run.json' if args.section6_only else 'full_run.json'
    for name,command in jobs:
        print('Running '+name,flush=True);start=time.monotonic()
        log=out/(name+'.log')
        with log.open('w') as stream:
            try:
                result=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=args.timeout)
                status='passed' if result.returncode==0 else 'failed'
                code=result.returncode
            except subprocess.TimeoutExpired:
                status,code='timeout',None
        item={'name':name,'status':status,'returncode':code,'seconds':round(time.monotonic()-start,2),
              'log':str(log.relative_to(ROOT))}
        report['jobs'].append(item)
        report['status']='running' if status=='passed' else 'failed'
        (out/report_name).write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(item),flush=True)
        if status!='passed':
            print(log.read_text()[-5000:]);raise SystemExit(1)
    report['status']='passed'
    report['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (out/report_name).write_text(json.dumps(report,indent=2)+'\n')
    print('All selected checks passed.',flush=True)

if __name__=='__main__':
    main()
