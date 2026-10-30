from run_v3_2 import *
import importlib.metadata as md,subprocess
rows=[]
for sec in M.index:
 p=O/'by_dataset'/M.loc[sec,'dataset_id']/sec;d=pd.read_parquet(p/'spot_regions_v3_2.parquet',columns=['shared_only_threshold_crossing','lesion_location','domain_component_id','domain_local_evidence_score','domain_shared_evidence_score','assignment_basis'])
 for did,z in d[d.shared_only_threshold_crossing].groupby('domain_component_id'):rows.append({'section_id':sec,'domain_component_id':did,'n_spots':len(z),'local_score':z.domain_local_evidence_score.iloc[0],'shared_score':z.domain_shared_evidence_score.iloc[0],'v3_2_locations':';'.join(sorted(z.lesion_location.unique())),'basis':';'.join(sorted(z.assignment_basis.unique()))})
 for file in [p/'small_roi_decisions.tsv']:
  a=pd.read_csv(file,sep='\t');a.rename(columns={'old_outer_segments':'old_roi_total_edge_segments'}).to_csv(file,sep='\t',index=False)
pd.DataFrame(rows).to_csv(O/'tables/shared_only_domain_crossings.tsv',sep='\t',index=False)
p=O/'tables/small_roi_decisions.tsv';a=pd.read_csv(p,sep='\t');a.rename(columns={'old_outer_segments':'old_roi_total_edge_segments'}).to_csv(p,sep='\t',index=False)
versions={x:md.version(x) for x in ['numpy','pandas','scipy','shapely','matplotlib','Pillow','pyarrow']};versions['python']=sys.version;versions['seed']=CFG['seed'];git=subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True);versions['repository_head']=git.stdout.strip() or None;versions['repository_status']='working project is not a git checkout; executed source files included' if git.returncode else 'git checkout';versions['author_SpatialBC_tree_commit']='f854197edb9730a9fc215ce1022a44474446f99c';(O/'software_versions.json').write_text(json.dumps(versions,indent=2))
review=[
 ('GSE319536_BC18','D3/D13 and lower/side unresolved tissue','844 remaining extent-unresolved spots; S1 Others and Tumor Stroma cannot be assigned as exact barcode truth'),
 ('GSE285715_T3','former outside plus weak local islands','1445 extent-unresolved; broad EPI and tumour source do not resolve macro extent; no Tumor-419a alias bridge'),
 ('GSE238145_B16','old L4 (5 spots), L6 (1 spot), former outside','small local-supported candidates retained; original outside lacks positive macro-external evidence'),
 ('GSE319536_BC01','author Normal Stroma vs project unresolved domains','published normal compartment visible qualitatively; absent per-spot mapping prevents precise transfer'),
 ('GSE319536_BC13','author Normal Stroma/Muscle/Necrosis compartments','do not equate all low malignant domains with outside; published image qualitative only'),
 ('Inhousedata_T4465935','former 4 candidate outer edges around v3.2 L3','QC missing channel fails external-connection check; now internal non-lesion compartment contact, no macro band'),
 ('Inhousedata_T4471438','old L2','1 weak shared-only spot withdrawn; other small local-supported candidate retained'),
 ('NPJ_sample1','old L8 and author sample mapping','1 weak shared-only spot withdrawn; paper alias bridge unavailable')]
pd.DataFrame(review,columns=['section_id','local_region','remaining_issue_or_resolution']).to_csv(O/'tables/local_cases.tsv',sep='\t',index=False)
