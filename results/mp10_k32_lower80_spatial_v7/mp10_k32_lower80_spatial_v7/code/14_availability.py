from common import *
import anndata as ad
mp=pd.read_csv(P/'mp10_k32_validation_v4/00_audit/top50_signatures.tsv',sep='\t');g0=set(pd.read_csv(R/'00_manifest/gene_manifest.tsv',sep='\t').gene);top20=set(pd.read_csv(R/'00_manifest/FROZEN_MP10_TOP20.tsv',sep='\t').gene);panel=['SREBF2','SREBF1','HSD17B7','FASN','ELOVL5','ADIPOR2','HMGCS1','ABHD3','INSIG1','ACSL3','ACSS2','EPCAM','KRT8','KRT18'];lookup=mp.groupby('gene').MP.agg(lambda x:';'.join(sorted(set(x))));genes=sorted(set(mp.gene)|set(panel));rows=[]
for sec in config()['sections']:
 a=ad.read_h5ad(P/'cell2location_preflight_v1/02_visium/anndata'/(sec+'.h5ad'),backed='r');available=set(a.var_names)
 for g in genes:rows.append({'section':sec,'dataset':sec.split('_')[0],'gene':g,'measured_feature_present':g in available,'in_G0':g in g0,'frozen_MP10_Top20':g in top20,'fixed_support_or_epithelial_panel':g in panel,'frozen_14MP_Top50_memberships':lookup.get(g,''),'rule':'exact feature ID; availability does not require nonzero expression in every spot'})
 a.file.close()
wt(pd.DataFrame(rows),R/'00_manifest/selected_section_gene_availability.tsv');say('FEATURE AVAILABILITY',len(rows))
