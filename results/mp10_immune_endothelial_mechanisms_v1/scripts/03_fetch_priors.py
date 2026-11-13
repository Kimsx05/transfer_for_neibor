import requests,pathlib,json,hashlib,time
P=pathlib.Path('resources/mp10_regulatory_priors');log=[]
def get(name,url):
 p=P/name
 try:
  if not p.exists():
   r=requests.get(url,timeout=(20,90),stream=True);r.raise_for_status()
   with open(str(p)+'.part','wb') as f:
    for z in r.iter_content(1024*1024):f.write(z)
   pathlib.Path(str(p)+'.part').rename(p)
  h=hashlib.sha256()
  with p.open('rb') as f:
   for z in iter(lambda:f.read(1024*1024),b''):h.update(z)
  x=dict(file=name,url=url,bytes=p.stat().st_size,sha256=h.hexdigest(),status='completed')
 except Exception as e:x=dict(file=name,url=url,status='technically_blocked',reason=str(e))
 log.append(x);print(x,flush=True);(P/'download_manifest.json').write_text(json.dumps(log,indent=2))
for n in ['lr_network_human_21122021.rds','ligand_target_matrix_nsga2r_final.rds']:
 get(n,'https://zenodo.org/records/7074291/files/'+n+'?download=1')
get('collectri_signed.csv','https://raw.githubusercontent.com/saezlab/CollecTRI/main/resources/CollecTRI_signed.csv')
get('collectri_omnipath.tsv','https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=references,sources&format=tsv')
get('progeny_model_human_full.rda','https://raw.githubusercontent.com/saezlab/progeny/master/data/model_human_full.rda')
get('dorothea_hs.rda','https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.rda')
get('human_tfs.txt','https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt')
get('motifs-v9-nr.hgnc-m0.001-o0.0.tbl','https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl')
get('hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather','https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather')
