import csv,re,hashlib
from pathlib import Path
from geo import D,I,put
T=D/'tables'
def rd(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def annotate():
 origins={'MP10':('neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt','original Top20 v5 calibration','metabolic transcriptional state'),'mTORC1':('mp10_hypothesis_followup_v2/tables/hallmark_gene_sets.tsv.gz','MSigDB 2026.1.Hs HALLMARK_MTORC1_SIGNALING','pathway activity proxy'),'SREBF1':('mp10_hypothesis_followup_v2/tables/scenic_corrected_regulon_edges.tsv','frozen v2 SREBF1(+)','regulon activity proxy'),'SREBF2':('mp10_hypothesis_followup_v2/tables/scenic_corrected_regulon_edges.tsv','frozen v2 SREBF2(+)','regulon activity proxy'),'Treg':('highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1/models/signatures.json','Treg_main 11 genes','bulk Treg composition proxy'),'Tlineage':('highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1/models/signatures.json','T_lineage 5 genes','bulk T-lineage composition proxy')}
 p=T/'variable_dictionary.tsv';var=rd(p)
 for z in var:
  name=z['variable']
  if name in origins:z['source'],z['version'],z['kind']=origins[name]
  elif '__' in name:
   n,b=name.rsplit('_',1);x,y=n.split('__');k=x if b=='A' else y;z['source'],z['version'],z['kind']=origins[k];z['version']+='; remove shared members of '+x+' and '+y
  else:z['version']='individual platform SOFT / GPL570 annotation / TCGA RSEM 2018';z['kind']='single gene RNA'
 put(p,var)
 for a,path,col in [('GSE32894',I/'GPL6947_extracted.soft','Symbol'),('GSE31684',D.parents[1]/'results/mp10_resistance_validation_v2_expanded/00_plan_audit/GPL570_annotation.tsv','Gene Symbol'),('GSE13507',I/'GPL6102.soft','Symbol')]:
  genes=set()
  with path.open(errors='replace') as f:
   for line in f:
    h=line.strip().split('\t')
    if h[0]=='ID' and col in h:break
   j=h.index(col)
   for line in f:
    if line.startswith('!platform_table_end'):break
    r=line.rstrip('\n').split('\t')
    if len(r)<=j:continue
    gs={re.sub(r'^PVRL([1-4])$',r'NECTIN\1',v.strip()) for v in re.split(r' */// *| *// *| *; *',r[j])}-{ '', '---','NA','null'}
    if len(gs)==1:genes.update(gs)
  p=T/'bulk'/a/'score_coverage.tsv';rows=rd(p)
  for z in rows:z['platform_annotated_n']=len(set(z['members'].split(';')) & genes)
  put(p,rows)
if __name__=='__main__':annotate()

