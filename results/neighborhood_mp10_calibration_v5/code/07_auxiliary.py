from common import *
import anndata as ad
from scipy.io import mmwrite

def main(out):
 P=out/'03_attribution_calibration';dm=pd.read_csv(P/'simulation_donors.tsv.gz',sep='\t');r=ad.read_h5ad(ROOT/'subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad',backed='r');refs=set(r.obs_names);source=set(r.obs['orig.ident'].astype(str));r.file.close();rows=[]
 for donor,g in dm.groupby('donor'):
  ids=set(';'.join(g.cell_ids).split(';'));hit=len(ids&refs);rows.append(dict(donor=donor,n_simulation_cells=len(ids),n_exact_reference_cell_matches=hit,donor_name_in_reference_origident=donor in source,reference_status='IN_REFERENCE' if hit or donor in source else 'DONOR_OVERLAP_UNVERIFIED',basis_status='original projection input source; not independent of MP basis',held_out=False))
 write_tsv(out,P/'donor_reference_overlap.tsv',rows)
 sec='CanDis_2026_Sample1';folder=V4/'05_source_validation/counts_cache'/sec;axis=json.loads((folder/'axis_manifest.json').read_text());Y=sparse.load_npz(folder/'composition_observed.npz').toarray()[:12];keys=axis['composition_ids'][:12];cm=pd.read_csv(V4/'05_source_validation/composition_only_pseudospots.tsv',sep='\t').set_index('bag_id').loc[keys];d=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet').set_index('spot_uid');a,C,s=params(sec);ids=a.obs_names.get_indexer(d.loc[cm.blueprint_uid,'barcode']);E=allocate(Y,C,np.asarray(s['w_sf'])[ids],s['m_g'],np.asarray(s['s_g_gene_add'])[a.obs['_scvi_batch'].to_numpy(int)[ids]]);allx=np.vstack([Y,E]);meta=[]
 for j,x in enumerate(allx):
  sc,_,_=auc_parts(x[None,:],axis['genes']);vals=x[x>0];meta.append(dict(cell_id=f'calibration_{j}',Top20_score=sc[0],input_type='integer_synthetic_observed' if j<len(Y) else 'fractional_estimated_epi',zero_fraction=(x==0).mean(),positive_tie_fraction=1-len(np.unique(vals))/len(vals) if len(vals) else np.nan))
 f=out/'01_score_fidelity';mmwrite(safe(out,f/'official_synthetic_input.mtx'),sparse.csc_matrix(allx.T));write_tsv(out,f/'official_synthetic_axes.tsv',meta);safe(out,f/'official_synthetic_genes.txt').write_text('\n'.join(axis['genes'])+'\n')
 write_json(out,out/'00_audit_and_spec/software_versions.json',{k:__import__('importlib.metadata',fromlist=['version']).version(k) for k in ['numpy','pandas','scipy','anndata','cell2location','scvi-tools','torch']});log('donor overlap and fractional comparison prepared')
if __name__=='__main__':main(OUT)
