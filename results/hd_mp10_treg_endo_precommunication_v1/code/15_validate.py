from common import *
from scipy.spatial import cKDTree
from pathlib import Path
import importlib.util
checks=[]
def ck(name,value,detail=''):
 checks.append(dict(check=name,pass_check=bool(value),detail=detail));assert value,(name,detail)
# One targeted artificial geometry test: center excluded, Epi included, Unknown outside conditional denominator, no component mixing.
xy=np.array([[0.,0.],[8,0],[16,0],[50,0],[51,0],[0,0]]);typ=np.array(['Epi','Epi','Unknown','Endo','T','T']);comp=np.array([1,1,1,1,1,2]);ix=np.where((np.linalg.norm(xy-xy[0],axis=1)<=50)&(comp==1)&(np.arange(6)!=0))[0];known=typ[ix][typ[ix]!='Unknown'];ck('physical_radius_center_unknown_component',list(ix)==[1,2,3] and set(known)=={'Epi','Endo'} and len(known)==2)
g=pd.read_csv(R/'tables/bin16_gene_conservation.tsv',sep='\t');ck('bin16_per_gene_conservation',np.array_equal(g.input_UMI,g.aggregated_UMI));ck('bin16_global_conservation',g.input_UMI.sum()==64528669)
d=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');ck('cell_IDs_unique',d.cell_id.is_unique);ck('original_analysis_QC_count',len(d)==93593);ck('original_analysis_UMI',d.nCount_RNA.sum()==60670091);ck('Top20_genes_frozen',len(TOP)==20 and set(TOP)==set('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()))
f=pd.read_csv(R/'tables/reference_model_features.tsv',sep='\t');ck('TIGIT_not_classifier_feature','TIGIT'not in set(f.gene));z=pd.read_csv(R/'tables/official_AUCell_validation.tsv',sep='\t');ck('official_AUCell_same_rank_formula',z.formula_equal.all() and np.max(np.abs(z.official_auc-z.same_rank_formula))<1e-12)
for r in [50,200]:
 n=pd.read_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{r}.parquet');ck(f'r{r}_all_Epi_centers',set(n.cell_id)==set(d.loc[d.epi_anchor,'cell_id']));pc=[c for c in n if c.startswith('P_')];cc=['count_'+c[2:]for c in pc];ck(f'r{r}_Epi_in_composition','P_Epithelial'in pc);ck(f'r{r}_count_denominator',np.array_equal(n[cc].sum(axis=1),n.n_known));ck(f'r{r}_conditional_sum',np.allclose(n.loc[n.n_known>0,pc].sum(axis=1),1));ck(f'r{r}_zero_known_not_imputed',(n.loc[n.n_known==0,'niche']==0).all());ck(f'r{r}_counts_not_UMI',n.n_neighbors.max()<10000)
 rng=np.random.default_rng(123);sample=n.iloc[rng.choice(len(n),min(12,len(n)),False)];types=[c[2:]for c in pc]
 for row in sample.itertuples():
  q=d[d.component==row.component];dist=np.hypot(q.x_um-row.x_um,q.y_um-row.y_um);q=q[(dist<=r+1e-7)&q.cell_id.ne(row.cell_id)];ct=q.type_supported.value_counts();expected=np.array([ct.get(t,0)for t in types]);actual=n.loc[n.cell_id.eq(row.cell_id),cc].to_numpy()[0];ck(f'r{r}_direct_geometry_{row.cell_id}',np.array_equal(expected,actual))
reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');mem=pd.read_csv(R/'tables/epithelial_region_membership.tsv.gz',sep='\t');ck('Epi_region_members_unique',mem.cell_id.is_unique);ck('Epi_region_no_TME',set(mem.cell_id)<=set(d.loc[d.mp10_epithelial_eligible,'cell_id']));ck('region_mapping_count',len(mem)==reg.n_epi.sum());ck('region_counts_conservation',reg.nCount_RNA.sum()==d.set_index('cell_id').loc[mem.cell_id,'nCount_RNA'].sum());ck('not_replicated_region_scores',reg.region_id.is_unique)
a=pd.read_csv(R/'tables/association_results.tsv',sep='\t');ck('Treg_NE_not_zero_effect',a.loc[a.predictor.str.startswith('Treg')|a.predictor.str.startswith('reverse'),'P'].isna().all());ck('freeze_precedes_association',(R/'niche_frozen.json').stat().st_mtime<(R/'models/association_execution.json').stat().st_mtime)
sp=importlib.util.spec_from_file_location('loader',R/'code/12_loaders.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);o=m.load_precommunication();ck('Python_loader_ID_alignment',o.index.is_unique and set(o.index)==set(d.cell_id));put(pd.DataFrame(checks),'validation_checks.tsv');print('VALIDATED',len(checks))
