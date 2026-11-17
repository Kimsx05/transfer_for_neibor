from common import *
import sys
s=sys.argv[1];gate();(R/s).mkdir(exist_ok=True)
if (R/s/'prepared.json').exists() and os.environ.get('FORCE')!='1':print('cached',s);raise SystemExit
x,d,ge,sy,path=load8(s)
# No pre-aggregation count filter. Only originally confirmed outside-tissue positions are excluded.
keep=d.in_tissue_original.isna().to_numpy() | (d.in_tissue_original.to_numpy()==1)
raw=x[:,np.isin(sy,PANELS['Epi'])];epi_n=np.asarray((raw>0).sum(1)).ravel()
anyepi=(d.type_blca.eq('Epithelial')|d.type_hpca.eq('Epithelial')).to_numpy()
epi=keep&(epi_n>=2)&(anyepi|(epi_n>=3))
# Preserve all supplied nodes and raw selected-gene evidence for plots and exact geometry.
cols=['bin_id','array_row','array_col','x_um','y_um','component','in_tissue_original','analysis_qc','nCount_RNA','nFeature_RNA','type_blca','type_hpca','reference_conflict']
b=d[cols].copy();b['eligible_pool']=keep;b['epi_support']=epi
for g in sorted(set(sum(PANELS.values(),[])+COMM+TOP)):
 b['raw_'+g]=np.asarray(x[:,sy==g].sum(1)).ravel()
b.to_parquet(R/s/'bin_evidence.parquet',index=False)
audit=[]
for kind,side,mask in [('T',32,keep),('T',16,keep),('Epi',64,epi)]:
 audit.append(summarize(x,d,sy,kind,side,mask,s));print(audit[-1],flush=True)
put(pd.DataFrame({'gene_id':ge,'symbol':sy}),s+'_gene_mapping.tsv')
info={'section':s,'counts_path':str(path),'counts_bytes':path.stat().st_size,'counts_mtime_ns':path.stat().st_mtime_ns,'input_bins':len(d),'input_genes':len(ge),'input_UMI':int(x.sum()),'original_tissue_mask':s in SECTIONS[:3],'coordinates_source':str(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet'),'audit':audit}
(R/s/'prepared.json').write_text(json.dumps(info,indent=2));print('PREPARED',s,flush=True)
