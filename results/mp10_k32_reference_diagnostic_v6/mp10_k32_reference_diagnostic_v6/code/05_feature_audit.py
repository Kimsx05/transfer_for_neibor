from common import *
O=R/'01_feature_audit';g=pd.read_csv(S/'02_inputs/gene_manifest.tsv',sep='\t').set_index('gene');inv=pd.read_csv(S/'02_inputs/spatial_actual_feature_inventory.tsv.gz',sep='\t');fm=pd.read_csv(P/'cell2location_full_run/v2_30type/01_manifest/01_spatial_manifest_frozen.csv');sections=inv[['dataset','section']].drop_duplicates();assert len(sections)==63
src=Path('/home/data/t070721/Metabolistic_scRNA/Multi2/4.3 NMF allgene/output_full_cc_seed5_v2');basis=pd.read_csv(src/'07_EpiMP_basis_program_weighted.csv').set_index('gene');top=pd.read_csv(src/'07_EpiMP_top50_genes.csv');markers=pd.read_csv(src/'07_EpiMP_marker_genes.csv');wt(markers,O/'frozen_marker_ranks.tsv');wt(top,O/'frozen_top50.tsv')
# Exact original feature-to-canonical relationship inherited from the frozen gene audit.
canonical={}
for gene,row in g.iterrows():
 for x in set(str(row.source_features).split(';'))|{gene}:canonical.setdefault(x,set()).add(gene)
def mapped(genes):
 out=set();unresolved=[]
 for x in genes:
  hit=canonical.get(x,set())
  if len(hit)==1:out.update(hit)
  else:unresolved.append(x)
 return out,unresolved
sets={}
for mp in basis.columns:
 sets[(mp,'full_positive_basis')]=set(basis.index[basis[mp]>0]);sets[(mp,'Top50')]=set(top.loc[top.EpiMP==mp,'gene']);sets[(mp,'Top20')]=set(top.loc[(top.EpiMP==mp)&(top['rank']<=20),'gene'])
keys=['HSD17B7','HMGCS1','KRT18','KRT8','EPCAM','SREBF2','FASN','ACSS2','INSIG1','ELOVL5'];rows=[];coverage=[];panels=[];g0=set(g.index[g.selected]);eligible=set(g.index[g.passes_TRAIN_filter & g.resolved_identity & ~g.is_MT])
feature_by={s:set(d.gene) for s,d in inv.groupby('section')}
for r in sections.itertuples():
 features=feature_by[r.section];frow=fm.loc[fm.section_id==r.section];hpath=frow.input_path.iloc[0] if len(frow) else 'UNRESOLVED'
 for gene in keys:
  present=gene in features
  rows.append(dict(dataset=r.dataset,section=r.section,assay='Visium-like image class; chemistry UNRESOLVED',panel_version='UNRESOLVED',gene=gene,original_feature_probe_ID=pd.NA,canonical_ID=g.loc[gene,'ensembl_gene_id'] if gene in g.index else pd.NA,original_measurable='UNKNOWN',processing_record='Preflight copies Spatial/counts gene rownames without gene filtering; earlier import/QC unavailable',existing_h5ad_retained=present,G0_retained=gene in g0,TRAIN_filter_pass=bool(g.loc[gene,'passes_TRAIN_filter']) if gene in g.index else False,evidence_file=str(hpath),script_evidence=str(P/'cell2location_preflight_v1/code/preflight/04_prepare_visium_inputs.R'),conclusion='PRESENT_IN_DERIVED_MATRIX_EXCLUDED_BY_63_SECTION_INTERSECTION' if present and gene not in g0 else ('PRESENT_IN_G0' if gene in g0 else 'RAW_ABSENCE_VS_PROCESSING_LOSS_UNRESOLVED'),status='PARTIAL'))
 for (mp,definition),original in sets.items():
  genes,unresolved=mapped(original);coverage.append(dict(dataset=r.dataset,assay_evidence='derived matrix; chemistry unresolved',section=r.section,MP=mp,definition=definition,n_frozen_original=len(original),n_unique_canonical=len(genes),n_unresolved_original=len(unresolved),n_matrix_present=len(genes&features),n_G0=len(genes&g0),fraction_frozen_exactly_mapped_and_present=len(genes&features)/len(original),unresolved_original=';'.join(sorted(unresolved))))
for ds,d in sections.groupby('dataset'):
 fs=set.intersection(*(feature_by[s] for s in d.section));genes=sorted(eligible&fs);pid='DERIVED_'+str(ds).replace('/','_');wt(pd.DataFrame({'gene':genes}),O/(pid+'_genes.tsv'));panels.append(dict(panel=pid,dataset=ds,sections=';'.join(d.section),n_sections=len(d),n_genes=len(genes),n_added_vs_G0=len(set(genes)-g0),added_key_genes=';'.join(x for x in keys if x in set(genes)-g0),status='PROVISIONAL_DERIVED_MATRIX_PANEL_REQUIRES_ASSAY_CONFIRMATION',basis='TRAIN-only permissive filter and exact within-dataset observed feature intersection; independent of predictions',gene_file=pid+'_genes.tsv',gene_hash=key(genes)))
wt(pd.DataFrame(rows),O/'gene_provenance_audit.tsv');wt(pd.DataFrame(coverage),O/'frozen_MP_gene_coverage_by_assay.tsv');wt(pd.DataFrame(panels),O/'G1_panel_candidates.tsv');wt(g.reset_index().groupby(['passes_TRAIN_filter','resolved_identity','is_MT','shared_all_eligible_sections','selected'],dropna=False).size().rename('n_genes').reset_index(),O/'G0_feature_loss_summary.tsv');wt(g.loc[keys].reset_index(),O/'key_gene_G0_loss.tsv')
js(R/'NEXT_G1_CONFIG.yaml',{'status':'NOT_RUN_BY_DESIGN','master_seed':123,'candidate_status':'PROVISIONAL; no assay-specific panel called validated without original chemistry/probe evidence','panels':panels,'before_training':['Confirm raw assay/panel version and actual feature/probe lineage for target sections','Freeze exact genes using MAIN TRAIN only plus confirmed assay measurability; remove MT and unresolved mappings','Matched S31 and S30 must use identical TRAIN members, gene order, priors and same DEV counts projected onto that panel','Retain existing frozen MP10 weights/Top20 and K32 membership; do not add substitutes','No selection using historical TEST or rescued-source/model performance'],'TEST_policy':'Already viewed; exploratory only, no new independent claim','planned_fits_this_run':0})
js(O/'FEATURE_AUDIT_STATUS.json',{'status':'PARTIAL','resolved':'Derived gene presence, G0 selection and preflight no-gene-filtering established','unresolved':'Raw chemistry/panel versions and pre-combined-object filtering lack original feature/probe evidence','searched_root':'/home/data/t070721/Metabolistic_scRNA','search_file':str(O/'raw_feature_search.txt'),'n_raw_annotation_files_found':len((O/'raw_feature_search.txt').read_text().splitlines()),'frozen_basis_file':str(src/'07_EpiMP_basis_program_weighted.csv'),'basis_sha256':sha(src/'07_EpiMP_basis_program_weighted.csv'),'Top20_definition':'rank<=20 from frozen 07_EpiMP_top50_genes.csv, not v4 differential-expression Top20','eligible_sections':63})
say('FEATURE_AUDIT_DONE',[(x,int(g.loc[x,'n_spatial_sections_present'])) for x in keys])
