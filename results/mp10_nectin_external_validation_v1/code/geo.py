from fetch import *
import re,collections
def read(p):
 md=collections.defaultdict(list)
 with gzip.open(p,'rt') as f:
  for r in csv.reader(f,delimiter='\t'):
   if not r:continue
   if r[0]=='!series_matrix_table_begin':break
   md[r[0]].append(r[1:])
 return md
def put(p,d):
 if not d:return
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,list(dict.fromkeys(k for r in d for k in r)),delimiter='\t');w.writeheader();w.writerows(d)
def run():
 audit=[]
 for a in accs:
  rows=[]
  for p in sorted(I.glob(a+'*series_matrix.txt.gz')):
   md=read(p)
   if a=='GSE32894':put(D/'config/extraction_32894.tsv',[dict(GSM=g,protocol=v) for g,v in zip(md['!Sample_geo_accession'][0],md['!Sample_extract_protocol_ch1'][0])])
   for j,g in enumerate(md['!Sample_geo_accession'][0]):
    def val(k):return md['!Sample_'+k][0][j]
    title=val('title');source=val('source_name_ch1');r=dict(accession=a,GSM=g,title=title,source=source,platform=val('platform_id'),organism=val('organism_ch1'))
    for ar in md['!Sample_characteristics_ch1']:
     if ':' in ar[j]:
      k,v=ar[j].split(':',1);r[k.lower().strip()]=v.strip()
    r['processing']=' | '.join(dict.fromkeys(ar[j] for ar in md['!Sample_data_processing'] if ar[j]))
    if a in accs[:3]:
     primary=a=='GSE13507' and source=='Primary bladder cancer'
     r['included']=int(a!='GSE13507' or primary)
     r['decision']='eligible tumor' if r['included'] else 'excluded: '+source
     r['stage_used']=r.get('tumor_stage',r.get('rc_stage',r.get('stage','')))
     st=r['stage_used'].upper();r['invasion']='MIBC' if re.search(r'T[234]',st) else 'NMIBC' if re.search(r'T(?:A|IS|1)',st) else 'unknown'
     r['T_category_invasion']=r['invasion']
     if a=='GSE13507' and r.get('invasiveness','').lower() in ('muscle invasive','non-muscle invasive'):
      r['invasion']='MIBC' if r['invasiveness'].lower()=='muscle invasive' else 'NMIBC'
     r['subtype_used']=r.get('molecular_subtype','')
     r['treatment_used']=r.get('prerc_chemo',r.get('systemic chemo','unknown'))
     r['patient_sample_id']=title.split('.')[0] if a=='GSE32894' else title.split('_')[0] if a=='GSE31684' else title.split()[-1].split('.')[0]
     r['primary_recurrent']='primary' if primary else 'unspecified tumor history' if a!='GSE13507' else source
    rows.append(r)
  if not rows:continue
  dest=D/'tables'/('bulk' if a in accs[:3] else 'perturbation')/a
  if a in accs[:3]:dest=dest/'cohort_manifest'
  put(dest/'sample_inclusion.tsv',rows)
  audit.extend([dict(accession=a,platform=r['platform'],source=r['source'],processing=r['processing']) for r in rows])
 put(D/'config/preprocessing.tsv',[dict(zip(audit[0],t)) for t in dict.fromkeys(tuple(r.values()) for r in audit)])
if __name__=='__main__':run()
