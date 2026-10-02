#!/usr/bin/env python3
from __future__ import annotations
import gzip,json,os,subprocess,sys,tempfile,time
from datetime import datetime,timezone
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
cases=[]
def run(name, args, reject=False):
    t=time.monotonic()
    try:
        p=subprocess.run([sys.executable,*args],cwd=ROOT,capture_output=True,text=True,
           timeout=45,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
        ok=(p.returncode!=0 and 'FAIL:' in p.stderr) if reject else p.returncode==0
        rec={'name':name,'command':['python',*args],'passed':ok,'returncode':p.returncode,
             'stdout':p.stdout,'stderr':p.stderr,'seconds':round(time.monotonic()-t,3)}
    except subprocess.TimeoutExpired:
        rec={'name':name,'passed':False,'reason':'45-second timeout'}
    cases.append(rec);print(('PASS ' if rec['passed'] else 'FAIL ')+name,flush=True)
for script in ['kazhdan_compressor_presentation','rigid_channel_builder','rigid_channel_instance']:
    run('baseline '+script,[f'verify/{script}.py'])
for variant in ['one_anchor','two_anchor']:
    run(variant,[f'verify/pruned_channel_check.py','--variant',variant])
for opt in [False,True]:
    for flag in ['coefficient','triangle','scope']:
        run('reject '+flag+(' optimized' if opt else ''),
            (['-O'] if opt else [])+['verify/pruned_channel_check.py','--perturb',flag],True)
for m in [2,10,100,1000]:
    re,im=Fraction(m*m-1,m*m+1),Fraction(2*m,m*m+1)
    ok=(re*re+im*im==1 and (1-re)**2+im*im==Fraction(4,m*m+1))
    cases.append({'name':f'rational phase identity m={m}','passed':ok,
                  'lazy_gap':str(Fraction(1,43*(m*m+1)))})
ok=(Fraction(13)*Fraction(25,12)<32 and 4096*19**2*44086==65188028416 and
    Fraction(1,447)**2 < Fraction(148534*41,(25*44086)**2))
cases.append({'name':'extraction constants','passed':ok})
with tempfile.TemporaryDirectory() as tmp:
    p=subprocess.run([sys.executable,'build_pruned.py','--emit',tmp],cwd=ROOT,
                      capture_output=True,text=True,timeout=45)
    ok=p.returncode==0; binding={}
    for name in ['one_anchor','two_anchor']:
        for tail in ['_channel.json.gz','_summary.json']:
            a=(Path(tmp)/(name+tail)).read_bytes();b=(ROOT/'data'/(name+tail)).read_bytes()
            if tail.endswith('.gz'):a,b=gzip.decompress(a),gzip.decompress(b)
            binding[name+tail]=a==b;ok=ok and a==b
    cases.append({'name':'deterministic re-emission','passed':ok,'bindings':binding})
passed=sum(c['passed'] for c in cases)
report={'executed_at_utc':datetime.now(timezone.utc).isoformat(),
        'counts':{'passed':passed,'total':len(cases)},'cases':cases,
        'scope':'Finite data and certificate arithmetic, not formal verification of the analytical proofs.',
        'global_gap_lower_bound':None,'external_rigidity_formalized':False}
(ROOT/'verification/quantitative_results.json').write_text(json.dumps(report,indent=2)+'\n', newline='\n')
print(f'RESULT {passed}/{len(cases)}',flush=True)
sys.exit(0 if passed==len(cases) else 1)
