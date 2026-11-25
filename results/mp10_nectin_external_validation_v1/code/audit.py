import sys,csv,hashlib,subprocess,platform,os,json,gzip
from pathlib import Path
from fetch import D,I,accs
from geo import put
from prepare import sha
T=D/'tables';root=D.parents[1];pr=root/'results/mp10_epithelial_treg_tcga_visualization_v1';v2=pr.with_name('mp10_epithelial_treg_tcga_visualization_v2')
def rd(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def rel(p):
 try:return str(p.relative_to(root))
 except:return str(p)
def run():
 rec=[];urls={}
 for r in csv.reader((D/'logs/download.tsv').open(),delimiter='\t'):
  if len(r)>=4 and r[3]!='0':urls[r[0]]=r[1]
 urls['GPL6947_extracted.soft']='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6947/soft/GPL6947_family.soft.gz'
 manifest=D/'config/resource_manifest.tsv'
 def reg(p,role,url=''):
  if p.is_file():rec.append(dict(path=rel(p),role=role,url=url,bytes=p.stat().st_size,sha256=sha(p),format='gzip' if p.suffix=='.gz' else p.suffix.lstrip('.')))
 for p in sorted(I.iterdir()):
  if p.suffix not in ['.part','.tmp'] and not p.name.startswith('mp10_'):reg(p,'official processed input',urls.get(p.name,''))
 for n,u in [('mp10_tcga_complete.txt','https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/blca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt'),('mp10_abs.tsv','https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5'),('mp10_clinical.txt','https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/blca_tcga_pan_can_atlas_2018/data_clinical_sample.txt')]:reg(Path('/tmp')/n if (Path('/tmp')/n).exists() else I/n,'TCGA cache',u)
 for p in ['neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt','highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1/models/signatures.json','results/mp10_resistance_validation_v2_expanded/00_plan_audit/GPL570_annotation.tsv','results/mp10_epithelial_treg_tcga_visualization_v1/scripts/build.R','results/mp10_epithelial_treg_tcga_visualization_v1/provenance/TCGA_sha256.txt','results/mp10_epithelial_treg_tcga_visualization_v2/plot_data/diagnostic_specimen_profiles.tsv']:reg(root/p,'frozen project input')
 for name in ['hallmark_gene_sets.tsv.gz','scenic_corrected_regulon_edges.tsv','sample_level_scores.tsv.gz','official_gene_sets.tsv.gz']:reg(root/'results/mp10_hypothesis_followup_v2/tables'/name,'frozen project input')
 for n in ['profiles.tsv','treg_logCPM.rds','perturbation_DE.tsv.gz','disjoint_control_sensitivity.tsv','original_CellChat.tsv']:reg(root/'results/mp10_treg_metabolic_mechanisms_v1/tables'/n,'frozen project input')
 for p in sorted((v2/'statistics').glob('TCGA_primary*.tsv')):reg(p,'original TCGA test')
 for mask in ['Full','M20']:reg(v2/'plot_data'/('TCGA_primary_program_'+mask+'_MP10_mTORC1.tsv'),'original TCGA scores')
 put(manifest,rec)
 orig=rd(T/'TCGA_original_Pearson.tsv');z=rd(T/'all_cohorts.tsv')
 qc=[]
 for s in z:
  if s['cohort']!='TCGA' or s['family_id']!='core':continue
  u=next(r for r in orig if r['source_file']=='TCGA_primary_program_Full_'+s['x']+'_'+s['y']+'.tsv' and r['mode']=='raw')
  delta=abs(float(s['r'])-float(u['r']));assert delta<1e-10;qc.append(dict(pair=s['x']+' / '+s['y'],original_r=u['r'],restored_r=s['r'],absolute_error=delta,original_q=u['q'],current_nine_pair_q=s['q']))
 put(D/'qc/TCGA_fidelity.tsv',qc)
 vars=[];tmp=rd(T/'bulk/TCGA/score_coverage.tsv')
 for x in tmp:
  vars.append(dict(variable=x['variable'],layer=x['mask'],kind='bulk program',input='TCGA log2(RSEM+1) / chip processed noncentered',source='frozen Top20 / Hallmark / scenic_corrected_regulon_edges / HD signatures.json',algorithm='GSVA ssGSEA alpha .25 normalize FALSE',members=x['members'],coverage='>=5 and >=50% of current set',provenance='reconstructed original' if x['variable'] in ['MP10','SREBF1','SREBF2','mTORC1'] else 'supplemental'))
 for x in ['NECTIN1','NECTIN2','NECTIN3','NECTIN4','PVR','TIGIT','PVRIG','CD226','CD3D','SREBF1_RNA','SREBF2_RNA','MTOR']:vars.append(dict(variable=x,layer='RNA',kind='bulk single gene',input='cohort-specific normalized expression',source='RSEM / GPL annotation',algorithm='unique annotated probes; per-sample median',members=x.replace('_RNA',''),coverage='measured finite pairwise',provenance='original NECTIN family or TF RNA; supplemental receptors and MTOR'))
 put(T/'variable_dictionary.tsv',vars)
 items=[('TCGA primary','mp10_epithelial_treg_tcga_visualization_v1/scripts/build.R + v2/statistics','原始40组配对版本及各自三类Pearson估计；核实407/397','原定义复用；外部3队列与跨队列FDR新增'),('TCGA family','mp10_nectin_family_rss_v1/statistics/correlations.tsv','已有Full/M20、Pearson/Spearman与纯度敏感性','新增外部复现'),('NECTIN2 / TIGIT / Treg','mp10_immune_endothelial_mechanisms_v1 + hypothesis_followup + HD signatures.json','已有通讯及空间代理。检索到的TCGA原表没有独立bulk Treg分数','同口径基因相关重建，11成员Treg明确为新增'),('paired Epi / Treg','mp10_hypothesis_followup_v1/v2','已有H05五项、NC/CY拆分及receiver G2M/OXPHOS','新增同分母单协变量竞争'),('Treg lipid and activation','mp10_treg_metabolic_mechanisms_v1/tables','已有SCAP/FASN/TCR/CD36/SREBF2方向投影与覆盖','复用。竞争模型仅重算去终点成员的控制读数'),('CellChat','mp10_epithelial_treg_tcga_visualization_v2 + mp10_treg_metabolic_mechanisms_v1','Treg→Epi RR 1.482；Epi→Treg RR 0.991；镜像配体受体','解释审计'),('真实干预','GSE236742 / GSE143803 / GSE267018','本轮实际完成相应RNA线别、培养条件和交互','新增'),('虚拟扰动','virtual_perturbation_manifest.json','只读接口，见 qc/virtual_interface.json','独立任务输出')]
 h=[dict(branch=a,source=b,previous=c,current=d) for a,b,c,d in items]
 put(T/'prior_analysis_audit.tsv',h);put(T/'reused_vs_new.tsv',h)
 for a in accs[:3]:
  rows=rd(T/'bulk'/a/'cohort_manifest/sample_inclusion.tsv');u=next(r for r in rows if r['included']=='1')
  put(T/'bulk'/a/'preprocessing_audit.tsv',[dict(accession=a,platform=u['platform'],original=u['processing'],program_input='reconstructed before centering' if a=='GSE32894' else 'official matrix',gene_input='official matrix',n_total=len(rows),n_included=sum(r['included']=='1' for r in rows))])
 if __name__=='__main__':print(len(rec),'input fingerprints')
if __name__=='__main__':
 run()
 from metadata import annotate
 annotate()

