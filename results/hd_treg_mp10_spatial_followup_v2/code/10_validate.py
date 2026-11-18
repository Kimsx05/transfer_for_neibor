from fu_common import *
sys.path.insert(0,str(R/"code"))
from load_result import load_followup,attach_region_evidence
import importlib.metadata as md,platform,subprocess
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)));assert ok,(name,detail)
for s in SECTIONS:
 d=load_followup(s);v=pd.read_parquet(V/s/'T32_regions.parquet').set_index('region_id');check(s+'_unique_ID',not d.region_id.duplicated().any());v=v.loc[d.region_id]
 check(s+'_score_QC_coordinates_unchanged',all(np.array_equal(d[c].to_numpy(),v[c].to_numpy(),equal_nan=True)for c in ['x_um','y_um','Treg_main','measurement_eligible']))
 check(s+'_no_label_change',d.old_evidence.equals(d.followup_evidence));check(s+'_no_communication_upgrade',not d.formal_communication_ready.any())
 out=attach_region_evidence(d[['region_id','x_um']],s);check(s+'_loader_aligned',out.fu2_followup_evidence.equals(d.followup_evidence))
 mm=load_followup(s,'membership');check(s+'_membership_unique',not mm.bin_id.duplicated().any()and mm.fu2_T32_region_id.isin(d.region_id).all())
 if s.startswith('MPBC'):
  b=pd.read_parquet(V/s/'bin_evidence.parquet');check(s+'_mask_unknown',b.in_tissue_original.isna().all() if 'in_tissue_original'in b else v.tissue_status.eq('unknown').all())
d=pd.read_parquet(R/'tables/Cervilla_paired_analysis.parquet');X=pd.read_parquet(R/'tables/design_r200_common.parquet');check('paired_design_ID_same_order',X.T_region_id.equals(d.T_region_id));Y=np.column_stack([zrank(d.MP10_bin),zrank(d.MP10_cell)]);b=np.linalg.lstsq(X.drop(columns='T_region_id'),Y,rcond=None)[0][-1];check('delta_equals_difference_model',np.isclose(np.diff(b)[0],np.linalg.lstsq(X.drop(columns='T_region_id'),Y[:,1]-Y[:,0],rcond=None)[0][-1]))
bs=pd.read_parquet(R/'tables/paired_bootstrap.parquet');check('all_paired_replicate_deltas_exact',np.allclose(bs.delta,bs.beta_cell-bs.beta_bin));check('both_source_depths_in_main_design',{'depth_bin','depth_cell'}<=set(X))
mm=pd.read_parquet(R/'tables/common_support_reconstruction_members.parquet');check('reconstructed_no_duplicate_source_per_T',not mm.duplicated(['T_region_id','source','member_id']).any());u=pd.read_csv(R/'tables/common_support_reconstructed_units.tsv',sep='\t');check('joint_score_eligibility',np.array_equal(u.paired_eligible,u.bin_eligible&u.cell_eligible));check('counts_integer_nonnegative',all((u[c]>=0).all()and(u[c]%1==0).all()for c in ['bin_nCount','cell_nCount']))
acct=pd.read_csv(R/'tables/common_support_count_accounting.tsv',sep='\t');check('common_support_per_gene_and_excluded_counts',len(acct)==494 and acct.per_gene_conserved.all()and np.array_equal(acct.old_source_raw_UMI,acct.geometry_retained_UMI+acct.geometry_excluded_UMI))
# Historical MP10 implementation and direct full-background average tied ranks agree on raw input.
x,nodes,genes,symbols,_=v1.load8('Cervilla');mem=pd.read_parquet(V/'Cervilla/Epi64_members.parquet');ev=pd.read_parquet(V/'Cervilla/Epi64_regions.parquet');idx=pd.Series(np.arange(len(nodes)),index=nodes.bin_id);sel=np.array([np.flatnonzero(symbols==g)[0]for g in v1.TOP]);diff=[]
for rr in ev[ev.measurement_eligible].head(4).itertuples():
 ids=mem.loc[mem.region_id.eq(rr.region_id),'bin_id'];v=np.asarray(x[idx.loc[ids]].sum(0)).ravel();rank=rankdata(-v,method='average');score=np.maximum(905-rank[sel],0).sum()/sum(905-np.arange(1,21));diff.append(abs(score-rr.MP10))
check('historical_MP1018085_full_rank_numerical',max(diff)<1e-12,max(diff));del x
parents=pd.read_parquet(R/'tables/mixture_parent_membership.parquet');check('reference_derived_profiles_partition_disjoint',parents.groupby('cell_id').partition.nunique().max()==1)
m=pd.read_parquet(R/'tables/expanded_mixtures.parquet');check('expanded_mixtures_count',len(m)==21600);check('realized_RNA_origin_mass',((m.realized_Treg_RNA_UMI<=m.realized_T_RNA_UMI)&(m.realized_T_RNA_UMI<=m.nCount)).all());check('native_v1_mixtures_exact_replay',pd.read_csv(R/'tables/mixture_replay_validation.tsv',sep='\t').passed.all())
r0=pd.read_csv(V/'tables/roi_manifest.tsv',sep='\t');r1=pd.read_csv(R/'tables/roi_manifest_display512.tsv',sep='\t');check('ROI24_identity_and_original_geometry',len(r1)==24 and r0[['section','roi','region_id','x0','x1','y0','y1']].equals(r1[['section','roi','region_id','x0','x1','y0','y1']]));check('display512_only',np.allclose(r1.display_x1-r1.display_x0,512,atol=1e-8,rtol=0)and np.allclose(r1.display_y1-r1.display_y0,512,atol=1e-8,rtol=0));check('all24figures_exist',all((R/'figures'/f'{r.section}_{r.roi}_512atlas.png').exists()for r in r1.itertuples()))
# Actual R loader, no Seurat object or matrix copy.
script=f"source({str(R/'code/load_result.R')!r}); d<-load_followup({str(R)!r},'Cervilla'); stopifnot(nrow(d)==19271,!anyDuplicated(d$region_id)); cat('R_LOADER_OK\\n')"
p=subprocess.run(['Rscript','--vanilla','-e',script],capture_output=True,text=True);check('R_loader_real_read',p.returncode==0,p.stdout+p.stderr)
put(pd.DataFrame(checks),'validation.tsv');versions=[]
for name in ['numpy','pandas','scipy','anndata','pyarrow','statsmodels','matplotlib','numba','PyYAML']:
 try:versions.append(dict(package=name,version=md.version(name)))
 except md.PackageNotFoundError:pass
versions.append(dict(package='Python',version=platform.python_version()));versions.append(dict(package='R',version=subprocess.run(['Rscript','--vanilla','-e','cat(as.character(getRversion()))'],capture_output=True,text=True).stdout));put(pd.DataFrame(versions),'software_versions.tsv')
print('VALIDATION_PASSED',len(checks),flush=True)
