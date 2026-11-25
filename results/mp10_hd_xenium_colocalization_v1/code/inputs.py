from base import *
import hashlib,importlib.metadata as md
hd=pd.read_csv(F/'tables/input_manifest.tsv',sep='\t');hd=hd[hd.section.isin(HD)]
items=[]
for row in hd.itertuples():
 s=row.section
 for role,p in [('8um counts',Path(row.counts_path)),('8um nodes + components',Path(row.nodes_path)),('Epi64 scores, eligibility, H05',F/s/'Epi64_regions.parquet'),('Epi64 raw membership',F/s/'Epi64_members.parquet'),('complete 16um bin membership',F/s/'T16_members.parquet'),('16um physical coordinates',F/s/'T16_measure.parquet')]:items.append((s,role,p))
mp=P/'mp10_cellchat_high_definition_sensitivity_v1/tables'
for role,p in [('RCTD reference',A/'hd_mp10_treg_endo_precommunication_v1/models/fine_reference.rds'),('RCTD reference membership',A/'hd_mp10_treg_endo_precommunication_v1/tables/reference_membership.tsv.gz'),('final raw score and reference membership',P/'mp10_sc_c2l_final_annotation_v1/tables/analysis_cells.tsv.gz'),('frozen sample thresholds + high denominators',mp/'specimen_quantiles_coverage.tsv'),('archived high membership check',mp/'membership/C30__ALL__H05.tsv.gz'),('30/32 correspondence',P/'subtype/v2/reference_training_v2/11_final_spatial_reference/00_32_to_30_label_mapping.csv'),('unsmoothed major markers',A/'cell_annotation_v2_1/code/programs.py'),('subtype markers',A/'cell_annotation_v2/code/marker_rules.py')]:items.append(('shared',role,p))
for r in pd.read_csv(O/'shared/xenium_inputs.tsv',sep='\t').itertuples():items.append((r.sample,r.role,Path(r.path)))
items.append(('shared','explicit Epi_all and H05 frozen barcode sets',P/'results/mp10_immune_endothelial_mechanisms_v1/tables/pseudobulk_membership.tsv.gz'))
for p in (O/'shared/msigdb').glob('*.json'):items.append(('shared','MSigDB v2026.2.Hs official JSON 2026-11-25',p))
put(pd.DataFrame([(s,r,str(p),p.stat().st_size,p.stat().st_mtime_ns,hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<35e6 else None) for s,r,p in items],columns=['sample','role','path','bytes','mtime_ns','sha256_small_files']),O/'shared/input_manifest.tsv')
m=hd[['section','dataset','specimen','patient_id']].rename(columns={'section':'sample'});m['platform']='HD';m['scope']='first round'
extra=pd.DataFrame([dict(sample=s,dataset='GSE326226',specimen='GSM'+str(9625559+i),patient_id='unknown',platform='Xenium',scope='first round') for i,s in enumerate(XE)]+[dict(sample=s,dataset='GSE318916',specimen=None,patient_id='unknown',platform='HD',scope='registered only') for s in ['MPBC1','MPBC2']]);put(pd.concat([m,extra]),O/'shared/samples.tsv')
(O/'shared/versions.json').write_text(json.dumps({k:md.version(k) for k in ['scipy','numpy','squidpy','pandas','anndata','pyarrow','h5py','numba','matplotlib','scikit-learn']},indent=2))
# HGNC lookup distinguishes symbol aliases from absent targets. Canonical matches never silently replace program members.
hm=pd.read_csv(P/'subtype/v2/reference_training_v2/04_gene_harmonization/01_gene_mapping_manifest_v2.csv',low_memory=False)
h=hm[['input_symbol','mapping_class','canonical_symbol','ensembl_gene_id']].drop_duplicates('input_symbol').set_index('input_symbol')
for s in XE:
 d=O/'xenium'/s;panel=pd.read_csv(d/'panel.tsv.gz',sep='\t').set_index('name');z=panel.join(h);z['reference_shared_exact']=z.index.isin((O/'shared/reference/genes.txt').read_text().splitlines());z['duplicated_name']=z.index.duplicated(keep=False);z['reference_canonical_collision']=z.canonical_symbol.duplicated(keep=False)&z.canonical_symbol.notna();put(z.reset_index(),d/'symbol_audit.tsv.gz')

