from common import *
import shutil
O=R/'02_inputs';old=V/'02_splits_and_genes';g=pd.read_csv(O/'gene_manifest.tsv',sep='\t',low_memory=False);sel=set(g.loc[g.selected,'gene']);mp=pd.read_csv(P/'mp10_k32_validation_v4/00_audit/top50_signatures.tsv',sep='\t')
wt(g[g.MP_memberships.notna()|g.fixed_support_panel],O/'MP14_MP10_fixed_panel_gene_audit.tsv')
a=pd.read_csv(old/'shared_gene_availability_by_dataset.tsv.gz',sep='\t');a=a.drop(columns=['provisional_common_gene_keep']);a['selected_scoped_MAIN_TRAIN']=a.gene.isin(sel);wt(a,O/'shared_gene_availability_by_dataset.tsv.gz');panel=set(g.loc[g.MP_memberships.notna()|g.fixed_support_panel,'gene']);wt(a[a.gene.isin(panel)],O/'MP14_MP10_fixed_panel_spatial_availability.tsv')
source_map={}
for z in g.itertuples():
 for x in set(str(z.source_features).split(';'))|{z.gene}:source_map.setdefault(x,[]).append(z.gene)
rows=[]
for z in mp.itertuples():
 targets=source_map.get(z.gene,[])
 rows.append({'MP':z.MP,'original_MP_input_gene':z.gene,'exact_source_mapped_genes':';'.join(targets),'selected_genes':';'.join(x for x in targets if x in sel),'retained':any(x in sel for x in targets),'mapping_rule':'exact original source feature or final canonical ID; no fuzzy matching'})
wt(pd.DataFrame(rows),O/'all_frozen_MP_input_feature_retention.tsv')
for f in ['patient_mapping_status.tsv','SPLIT_CONFIG.yaml']:shutil.copyfile(old/f,R/'01_scope'/f)
shutil.copyfile(old/'spatial_actual_feature_inventory.tsv.gz',O/'spatial_actual_feature_inventory.tsv.gz')
m=pd.read_csv(R/'01_scope/actual_cell_split_manifest.tsv.gz',sep='\t',low_memory=False)
wt(m.groupby(['split','new_reference_label'],observed=True).size().rename('n_cells').reset_index(),R/'01_scope/split_support_by_subtype.tsv')
wt(m.groupby(['split','dataset','new_reference_label'],observed=True).size().rename('n_cells').reset_index(),R/'01_scope/split_support_by_dataset.tsv')
print('Support panel:',g[g.fixed_support_panel][['gene','selected']].to_dict('records'))
print('MP source inputs retained:',sum(z['retained'] for z in rows),'/',len(rows))
# Identity audits are scope-independent; original files remain read-only.
for f in ['duplicate.tsv','annotation_conflict.tsv','unmatched_v4.tsv.gz','unmatched_old_epithelial.tsv.gz','old_compound_barcode_duplicates.tsv','v4_compound_barcode_duplicates.tsv','K32_retention.tsv']:
 shutil.copyfile(V/'01_label_crosswalk'/f,R/'01_scope'/('reused_identity_'+f))
scope=pd.read_csv(R/'01_scope/reference_scope_manifest.tsv.gz',sep='\t',low_memory=False)
wt(scope.groupby(['old_reference_label','K32_status'],dropna=False).size().rename('n_cells').reset_index(),R/'01_scope/old_subtype_by_K32_status.tsv')
provenance=m.groupby(['new_reference_label','new_major'],dropna=False).size().rename('n_cells').reset_index().rename(columns={'new_major':'original_source_major'})
wt(provenance,R/'01_scope/new_subtype_original_major_provenance.tsv')
provenance['reference_major']=np.where(provenance.new_reference_label.isin(['Epi_K32','Epi_other']),'Epithelial',provenance.original_source_major)
wt(provenance.groupby(['new_reference_label','reference_major']).n_cells.sum().reset_index(),R/'01_scope/new_subtype_major_mapping.tsv')
rows=[]
for name,mask in [('C14',m.frozen_cluster.eq('C14')),('C7',m.frozen_cluster.eq('C7')),('C17',m.frozen_cluster.eq('C17')),('UCell_high_nonK32',m.UCell_high.eq(True)),('NMF_high_nonK32',m.NMF_high.eq(True))]:
 t=m[mask&m.new_reference_label.eq('Epi_other')].groupby(['split','dataset','sample','source_group']).size().rename('n_cells').reset_index();t['stratum']=name;rows.append(t)
wt(pd.concat(rows,ignore_index=True),R/'01_scope/Epi_other_difficult_strata_support.tsv')
coverage=[];used_manifests=[]
for part in ['DEV','TEST']:
 mem=pd.read_csv(R/f'04_pseudospots/{part}/cell_membership_RNA.tsv.gz',sep='\t');seen=set(mem.cell_id);hold=m[m.split.eq(part)].copy();hold['used_in_any_pseudospot']=hold.cell_id.isin(seen);used_manifests.append(hold[['cell_id','split','dataset','sample','source_group','new_reference_label','used_in_any_pseudospot']])
 t=hold.groupby(['split','dataset','sample','source_group','new_reference_label']).agg(n_eligible=('cell_id','size'),n_unique_simulated=('used_in_any_pseudospot','sum')).reset_index();t['coverage']=t.n_unique_simulated/t.n_eligible;coverage.append(t)
wt(pd.concat(coverage,ignore_index=True),R/'01_scope/actual_heldout_cell_simulation_coverage.tsv')
print(pd.concat(coverage).query('new_reference_label=="Epi_K32"').groupby('split')[['n_eligible','n_unique_simulated']].sum().to_string())

wt(pd.concat(used_manifests,ignore_index=True),R/'01_scope/heldout_cell_simulation_use_manifest.tsv.gz')
status=[]
for part in ['DEV','TEST']:
 spots=pd.read_csv(R/f'04_pseudospots/{part}/manifest.tsv.gz',sep='\t')
 for dataset in sorted(m.dataset.unique()):
  z=spots[spots.dataset.eq(dataset)];k=m[m.split.eq(part)&m.dataset.eq(dataset)&m.new_reference_label.eq('Epi_K32')].groupby('source_group').size();pos=z[z.N_K32.gt(0)]
  state='NO_PSEUDOSPOTS_IN_SPLIT' if z.empty else ('NEGATIVE_ONLY_NO_POSITIVE_RECOVERY_EVALUATION' if pos.empty else 'POSITIVE_RECOVERY_EVALUABLE_WITH_LIMITED_UNVERIFIED_SOURCE_GROUPS')
  status.append({'split':part,'dataset':dataset,'n_pseudospots':len(z),'n_source_groups':z.source_group.nunique(),'n_positive_pseudospots':len(pos),'n_positive_source_groups':pos.source_group.nunique(),'n_eligible_K32_groups_ge30':int((k>=30).sum()),'status':state,'patient_status':'conservative source groups; not independently verified patients'})
wt(pd.DataFrame(status),R/'01_scope/evaluation_support_status.tsv')
