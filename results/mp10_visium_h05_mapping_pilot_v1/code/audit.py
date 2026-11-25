from common import *
import platform, scipy, pyarrow, matplotlib, h5py
records=read(O/'input_manifest.tsv').to_dict('records')
out=[]
for s,r in M.iterrows():
 d=read(O/f'spots/{s}.input.tsv.gz').set_index('barcode')
 p=Path(r.counts_path).parent/'spot_metadata.tsv'
 old=read(p).set_index('barcode')
 c=pd.read_csv(r.c2l_mean_path,index_col=0)
 for f in [old,c]:
  assert f.index.is_unique and set(f.index)==set(d.index)
 old=old.loc[d.index];c=c.loc[d.index]
 assert c.shape[1]==30 and set(old.section_id)=={s}
 frac=c.Epithelial/c.sum(axis=1,min_count=30)
 valid=(d.center_qc & d.lesion_location.ne('quality_excluded'))
 fields={k:np.nanmax(np.abs(old[k]-d[k])) for k in ['x','y','array_row','array_col','EpiMP10_Top20_AUCell']}
 assert all(v<1e-7 for v in fields.values())
 assert np.allclose(frac[valid],d.loc[valid,'Epithelial_fraction'],atol=1e-9)
 recipe=old.analysis_eligible_legacy & old.EpiMP10_Top20_AUCell.notna() & d.lesion_location.ne('quality_excluded')
 assert np.array_equal(recipe,d['query'])
 out.append(dict(section_id=s,legacy_query_recipe_matches=True,n=len(d),c2l_cols=len(c.columns),max_abs_fraction_difference=(frac[valid]-d.loc[valid,'Epithelial_fraction']).abs().max(),**fields))
 for role,path in [('independent_coordinate_score_audit',p),('mean_all30_audit',Path(r.c2l_mean_path))]:
  records.append(dict(role=role,section=s,path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
for role,path in [('H05_original_thresholds',R/'mp10_cellchat_high_definition_sensitivity_v1/tables/specimen_quantiles_coverage.tsv'),('confirmed_patient_pair',A.parent/'qc/patient_relationships.json'),('K32_historical_outcome',K/'REPORT.html'),('frozen_query_recipe',R/'mp10_paper_distance_full_v1/scripts/distance.py')]:
 records.append(dict(role=role,section='provenance',path=str(path),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
save(pd.DataFrame(records).drop_duplicates(['role','section','path']),'input_manifest.tsv')
save(out,'coordinate_score_audit.tsv')
(O/'versions_python.txt').write_text('\n'.join([sys.version,platform.platform()]+[f'{p.__name__} {p.__version__}' for p in [np,pd,scipy,pyarrow,matplotlib,h5py]])+'\n')

