from pathlib import Path
import sys,os,json,gzip,tarfile,io,re,hashlib
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');B=P/'results/mp10_resistance_validation_v1';S=B/'00_audit/sources'
sys.path=[p for p in sys.path if 'dist-packages' not in p];sys.path[:0]=[str(P/'spatial_ecology_top20_v2/runtime'),str(P/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
import pandas as pd,numpy as np
cfg=json.loads((B/'config/frozen_v1.json').read_text());assert hashlib.sha256((B/'config/frozen_v1.json').read_bytes()).hexdigest()==(B/'config/frozen_v1.sha256').read_text().strip()
def save(x,p):x.to_csv(B/p,sep='\t',index=False,na_rep='NA')
# Archive contents are individual processed gene counts, read directly without unsafe extraction.
t=tarfile.open(S/'GSE212810_RAW.tar');data=[]
for m in t.getmembers():
 assert m.isfile() and re.match(r'GSM\d+_Case_\d+.txt.gz$',m.name)
 d=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())),sep='\t');assert len(d.columns)==2
 data.append(d.set_index(d.columns[0]).rename(columns={d.columns[1]:m.name.split('_')[0]}))
assert all(x.index.equals(data[0].index) for x in data)
x=pd.concat(data,axis=1);x.index.name='original_id';x.to_csv(B/'03_clinical/GSE212810_counts_input.tsv.gz',sep='\t')
# Supplement table: first 9 rows are headers; explicit author column mapping.
r=pd.read_csv(S/'AMC_S6A_raw.tsv',sep='\t',header=None,skiprows=9)
cols={0:'patient_author',42:'case',43:'subtype',52:'response_table',53:'sex',58:'cT',62:'regimen_code',63:'cycles',64:'procedure',65:'ypT_code',66:'ypN_code'}
r=r[list(cols)].rename(columns=cols);r['case']=pd.to_numeric(r['case'].astype(str).str.extract(r'Case_#(\d+)')[0],errors='coerce');r=r[r.case.notna()].copy();r['case']=r.case.astype(int);assert r.case.is_unique
m=pd.read_csv(B/'00_audit/GSE212810_GSM_metadata.tsv',sep='\t');m['case']=m.title.str.extract(r'#(\d+)').astype(int)
m=m.merge(r,on='case',how='left',validate='one_to_one');m['author_good']=(m.patholresponse=='Yes').astype(int);m['regimen']=pd.to_numeric(m.regimen_code,errors='coerce').map({1:'GP_cisplatin',2:'GCb_carboplatin',3:'ddMVAC_cisplatin'})
m['pretreatment']=True;m['evidence']='GEO GSM + PMC10591055 Table S6A, case-number crosswalk footnote'
save(m,'03_clinical/GSE212810_clinical_input.tsv')
print('AMC regimen',m.regimen.value_counts(dropna=False).to_dict(),'matched',m.patient_author.notna().sum(),'label conflicts',((m.author_good!=pd.to_numeric(m.response_table,errors='coerce'))&m.response_table.notna()).sum())
# Cell-line mapping from GSM descriptions, no inference from expression values.
m=pd.read_csv(B/'00_audit/GSE171537_GSM_metadata.tsv',sep='\t');m['matrix_id']=m.description.str.split(' | ',regex=False).str[0];m['replicate']=m.title.str.extract(r'_([ABC])$');m['state']=m.resistance.map({'Cis':'CisR','Gem':'GemR','Gem_and_Cis':'GemCisR','Parental':'parental','parental':'parental'});m['cell_line']=m['cell line'];print(m[['resistance','state']].drop_duplicates().to_string(index=False));assert m.matrix_id.is_unique
m['replicate_unit']='RNA collected from 3 separate passages of same derived line; not independent resistance induction';m['washout']='several passages without antibiotics/gemcitabine/cisplatin (GSM growth protocol)';m['batch']='not separately documented';save(m,'02_cell_lines/sample_mapping.tsv')
# Local metadata audit only.
m=pd.read_csv(Path(cfg['inputs']['epi_cache'])/'metadata.tsv.gz',sep='\t',low_memory=False);f=pd.read_parquet(cfg['inputs']['epi_score']);print('Epi',m.shape,'fidelity',f.shape,'cycle',[c for c in m if 'Score' in c]);print('fidelity columns',list(f.columns));assert m.cell_id.is_unique and f.cell_id.is_unique
save(m.groupby(['dataset','Type'],dropna=False).size().reset_index(name='n_cells'),'00_audit/epi_dataset_source_counts.tsv')
p=pd.read_csv(cfg['inputs']['patient_map'],sep='\t');save(p,'00_audit/epi_sample_patient_mapping.tsv')
print('patient',p.patient_status.value_counts().to_dict());print('all preparation finished')
