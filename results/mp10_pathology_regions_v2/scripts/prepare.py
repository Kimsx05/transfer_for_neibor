import sys,json
from pathlib import Path
R=Path(__file__).resolve().parents[2];O=R/'mp10_pathology_regions_v2';E=R/'spatial_ecology_top20_v2';sys.path.insert(0,str(E/'runtime'))
import pandas as pd,numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
m=pd.read_csv(E/'00_inputs/section_manifest.tsv',sep='\t');im=pd.read_csv(E/'00_inputs/histology_manifest.tsv',sep='\t').set_index('section_id');d=pd.read_parquet(E/'00_inputs/frozen_spot_inputs.parquet');c=pd.read_parquet(E/'03_copykat_geography/spot_geography.parquet');idn=pd.read_csv(E/'patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv',sep='\t').set_index('section_id');pilot=['Inhousedata_T4413900','CanDis_2026_Sample8','GSE238145_B13','GSE285715_T3','GSE285715_N1','GSE171351_A1','GSE319536_BC01','NPJ_sample1','GSE238145_B10'];rows=[]
for sec,g in d.groupby('section_id',sort=False):
 g=g.reset_index(drop=True);ds=g.dataset.iloc[0];dest=O/'by_dataset'/ds/sec;dest.mkdir(parents=True,exist_ok=True);old=E/'03_copykat_geography/by_dataset'/ds/sec;counts=old/'input/raw_counts_gene_by_spot.mtx';counts=counts if counts.exists() else Path(str(counts)+'.gz');assert counts.exists()
 g=g[['section_id','barcode','x','y','array_row','array_col','epi_fraction_gate']];g['SpaCET_id']=g.array_row.astype(str)+'x'+g.array_col.astype(str);g['dataset_id']=ds;g['patient_id']=idn.loc[sec,'patient_uid'];g['source_context']='adjacent_normal_source' if sec.startswith('GSE285715_N') else 'tumor_study_specimen_not_region_diagnosis';g=g.merge(c[c.section_id==sec][['barcode','native_copykat_label','baseline_confidence','profile_available','source_adapted_CNA_cluster','spot_MSQ','in_automatic_baseline']],on='barcode',validate='one_to_one');g.to_csv(dest/'input_spots.tsv',sep='\t',index=False,na_rep='NA')
 lookup={(r.array_row,r.array_col):i for i,r in g.iterrows()};pairs=[]
 for (row,col),i in lookup.items():
  for dr,dc in [(0,2),(1,1),(1,-1)]:
   j=lookup.get((row+dr,col+dc))
   if j is not None:pairs.extend([(i,j),(j,i)])
 a=np.array(pairs);adj=coo_matrix((np.ones(len(a)),(a[:,0],a[:,1])),shape=(len(g),len(g))).tocsr();_,cc=connected_components(adj);two=adj+adj@adj;two.setdiag(0);two.eliminate_zeros();co=two.tocoo();edges=pd.DataFrame({'from':co.row+1,'to':co.col+1});edges.to_csv(dest/'banksy_two_ring_edges.tsv',sep='\t',index=False);g['tissue_component']=cc;g.to_csv(dest/'input_spots.tsv',sep='\t',index=False,na_rep='NA')
 rows.append(dict(dataset_id=ds,section_id=sec,n_spots=len(g),counts_path=str(counts),genes_path=str(old/'input/genes.txt'),barcodes_path=str(old/'input/barcodes.txt'),image_path=im.loc[sec,'image_path'],cna_path=str(old/'final_CNA_profile_original_barcodes.npz'),cluster_path=str(old/'source_adapted_clusters_original_barcodes.tsv'),c2l_mean_path=str(R/'cell2location_full_run/v2_30type/02_sections'/sec/'abundance_mean.csv.gz'),c2l_q05_path=str(R/'cell2location_full_run/v2_30type/02_sections'/sec/'abundance_q05.csv.gz'),output_dir=str(dest),pilot=sec in pilot,source_context=g.source_context.iloc[0]))
pd.DataFrame(rows).to_csv(O/'tables/section_manifest.tsv',sep='\t',index=False);(O/'tables/pilot_sections.json').write_text(json.dumps(pilot,indent=2));print(len(rows),'sections',sum(x['n_spots'] for x in rows),'spots; profiles existing',sum(Path(x['cna_path']).exists() for x in rows))
