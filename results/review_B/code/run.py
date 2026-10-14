from common import *
import argparse,subprocess

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['all','prepare','composition','association','figures','report'],default='all');ap.add_argument('--approval');ap.add_argument('--verify-only',action='store_true');args=ap.parse_args()
 if args.verify_only:
  subprocess.run(['/usr/bin/python3.12',str(OUT/'code/validate.py')],check=True)
  if (OUT/'reports/review_B.zip').exists():subprocess.run(['/usr/bin/python3.12',str(OUT/'code/validate.py'),'--delivery'],check=True)
  log('VERIFY_ONLY_PASS');return
 from prepare import prepare,rings
 prepare()
 if args.approval:
  supplied=pd.read_csv(args.approval,sep='\t');frozen=pd.read_csv(OUT/'00_approval/approved_gate_v7.tsv',sep='\t');cols=['candidate_id','scope','niche','celltype','accepted_label','decision','approved_scope'];assert supplied.candidate_id.is_unique and len(supplied)==1112;assert supplied.decision.eq('GO').all();assert supplied.sort_values('candidate_id')[cols].reset_index(drop=True).equals(frozen.sort_values('candidate_id')[cols].reset_index(drop=True)),'APPROVAL_SCOPE_CONFLICT'
 if args.stage in ['all','prepare']:rings()
 if args.stage in ['all','composition']:
  from composition import summarize,statistics
  summarize();statistics()
 if args.stage in ['all','association']:
  from association import section_associations,aggregate
  section_associations();aggregate()
 if args.stage in ['all','figures']:subprocess.run(['/usr/bin/python3.12',str(OUT/'code/launch_plots.py')],check=True)
 if args.stage in ['all','report']:
  subprocess.run(['/usr/bin/python3.12',str(OUT/'code/validate.py')],check=True)
  subprocess.run(['/usr/bin/python3.12',str(OUT/'code/report.py')],check=True)
  subprocess.run(['/usr/bin/python3.12',str(OUT/'code/validate.py'),'--delivery'],check=True)
  finalize()

def finalize():
 checkpoint('M06_DELIVERY_VERIFIED');s=pd.read_csv(OUT/'07_synthesis/scope_summary.tsv',sep='\t');v=pd.read_csv(OUT/'reports/validation_results.tsv',sep='\t');js('reports/execution_status.json',dict(status='PHASE_B_COMPLETE_AWAITING_RESULTS_REVIEW',n_pairs=139,n_scopes=8,n_approved_scope_pairs=1112,n_composition_figures=1112,n_MP10_figures=1344,n_structure_and_overview_figures=16,PNG_and_vector_PDF=True,validation_passed=int(v.status.eq('PASS').sum()),coverage=s.to_dict('records'),phase_B_real_computation=True,new_MP10_scoring=False,new_model_training=False,communication_analysis=False))
 # Refresh report/validation/status entries inside existing ZIP without duplicate archive names.
 import zipfile
 archive=OUT/'reports/review_B.zip';tmp=archive.with_suffix('.tmp.zip');replacement={str(p.relative_to(OUT)):p for p in (OUT/'reports').iterdir() if p.suffix in ['.md','.tsv','.json'] and p.name!='delivery_manifest.tsv'}
 for folder in ['00_audit','07_synthesis','code']:
  replacement.update({str(p.relative_to(OUT)):p for p in (OUT/folder).rglob('*') if p.is_file() and p.suffix!='.pyc' and '__pycache__' not in str(p)})
 with zipfile.ZipFile(archive) as src,zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as dst:
  for name in src.namelist():
   if name not in replacement and name!='reports/delivery_manifest.tsv':dst.writestr(name,src.read(name))
  for name,p in replacement.items():dst.write(p,name)
 tmp.replace(archive)
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None and len(z.namelist())==len(set(z.namelist()))
 tab('reports/delivery_manifest.tsv',[dict(path=p.name,size_bytes=p.stat().st_size,sha256=sha(p)) for p in sorted((OUT/'reports').glob('*.zip'))]);log('PHASE_B_COMPLETE_AWAITING_RESULTS_REVIEW')
if __name__=='__main__':main()
