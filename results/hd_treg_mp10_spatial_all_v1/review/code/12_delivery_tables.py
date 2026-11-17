from common import *
import sys,platform,importlib.metadata,re
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text())
info={
'Cervilla':('Zenodo17999961','visiumhd_Bladder','unknown','original_registered_HE','fullres pixels * validated mpp'),
'D1_Johnson':('GSE300692','GSM9066380/D1','unknown','original_registered_HE_crop','fullres pixels * validated mpp; retain crop origin'),
'UC_SCC':('GSE299573','GSM9040999/UC_SCC_20241126','unknown','downsampled_hires_only','fullres pixels * .3779763964 um/pixel'),
'MPBC1':('GSE318916','GSM9505917/mpbc4423','unknown','none','x=8*(col8+.5), y=8*(row8+.5); original 0-based array origin'),
'MPBC2':('GSE318916','GSM9505918/mpbc4538','unknown','none','x=8*(col8+.5), y=8*(row8+.5); original 0-based array origin')}
manifest=[];summary=[];roirows=[];geneinfo=[]
roi=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');nc=pd.read_csv(R/'tables/neighborhood_coverage.tsv',sep='\t')
for s in SECTIONS:
 identity=json.loads((R/s/'prepared.json').read_text());dataset,specimen,patient,he,coord=info[s]
 manifest.append(dict(section=s,dataset=dataset,specimen=specimen,patient_id=patient,patient_independence='not established across these HD sources',counts_path=identity['counts_path'],nodes_path=identity['coordinates_source'],input_bins=identity['input_bins'],genes=identity['input_genes'],input_UMI=identity['input_UMI'],resolution_um=8,coordinate_definition=coord,coordinates_unit='um',coordinate_direction='x right; y down; image orientation retained or relative grid',original_tissue_mask=s in SECTIONS[:3],HE=he,mask_limitation='input support is not histological tissue area'if s.startswith('MPBC')else 'supplied coverage smaller than possible image extent; unprovided area remains missing',reference='final30type frozen sample membership; old BLCA/HPCA background retained'))
 for branch in ['T32','T16','Epi64']:
  path=R/s/(branch+'_regions.parquet');d=pd.read_parquet(path)
  for name,genes in [('T_lineage',cfg['signatures']['T_lineage']),('Treg_main',cfg['signatures']['Treg_main']),('MP10',TOP),('Treg_core',['FOXP3','IL2RA','CTLA4','IKZF2'])]:
   cols=['raw_'+g for g in genes if 'raw_'+g in d]
   d[name+'_detected_genes']=[';'.join(g for g,c in zip([c[4:]for c in cols],vals)if c>0)for vals in d[cols].to_numpy()]
   d[name+'_raw_UMI']=d[cols].sum(axis=1)
  d['measurement_source']='supplied_8um_raw_sum';d['coordinate_status']='relative_array_tissue_unverified'if s.startswith('MPBC')else 'registered_um';d.to_parquet(path,index=False)
  if branch=='T32':t=d
  if branch=='Epi64':e=d
  geneinfo.append(dict(section=s,branch=branch,MP10_present=sum('raw_'+g in d for g in TOP),Treg_main_present=sum('raw_'+g in d for g in cfg['signatures']['Treg_main']),T_lineage_present=sum('raw_'+g in d for g in cfg['signatures']['T_lineage']),gene_background=18085))
 row=dict(section=s,T_measurement_regions=int(t.measurement_eligible.sum()),Treg_like_regions=int(t.evidence.eq('Treg_like_multigene').sum()),Epi_MP10_regions=int(e.measurement_eligible.sum()),credible_single_Treg_cells=0,mask_status='unknown'if s.startswith('MPBC')else 'original',patient_independence='unknown',status='completed; partially evaluable')
 for rad in [50,200]:
  v=nc[(nc.section==s)&nc.T_grid_um.eq(32)&nc.radius_um.eq(rad)&nc.source_rule.eq('exclude_shared')].iloc[0];row[f'comparable_T_regions_{rad}']=int(v.comparable_T_regions)
 summary.append(row)
 for r in roi[roi.section==s].itertuples():
  z=t.set_index('region_id').loc[r.region_id];er=e[e.x_um.between(r.x0,r.x1)&e.y_um.between(r.y0,r.y1)&e.measurement_eligible]
  roirows.append(dict(section=s,roi=r.roi,region_id=r.region_id,evidence=z.evidence,T_UMI=z.T_UMI,T_detected_genes=z.T_lineage_detected_genes,Treg_main_genes=z.Treg_main_detected_genes,Treg_core_genes=z.Treg_core_detected_genes,Treg_main=z.Treg_main,Tfh=z.Tfh,Activated_T=z.Activated_T,MP10_regions_in_ROI=len(er),interpretation='regional multigene evidence; no cell identity inferred; empty Epi panel means no eligible Epi score in displayed window'))
put(pd.DataFrame(manifest),'input_manifest.tsv');put(pd.DataFrame(summary),'section_summary.tsv');put(pd.DataFrame(roirows),'ROI_expression_evidence.tsv');put(pd.DataFrame(geneinfo),'gene_coverage.tsv')
versions={'Python':sys.version,'platform':platform.platform()}
for pkg in ['numpy','pandas','scipy','scikit-learn','anndata','h5py','numba','matplotlib','pyarrow','PyYAML']:
 try:versions[pkg]=importlib.metadata.version(pkg)
 except:versions[pkg]='not installed'
versions['UCell']='2.14.0 (native score core; R session in logs/UCell_sessionInfo.txt)';versions['SpottedPy']='0.1.7 source helper, not full package installation'
(R/'software_versions.json').write_text(json.dumps(versions,indent=2))
resources=[]
for f in (R/'logs').glob('*.log'):
 s=f.read_text(errors='replace');rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',s);wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',s)
 if rss or wall:resources.append(dict(step=f.stem,peak_RSS_KiB=int(rss.group(1))if rss else np.nan,wall=wall.group(1)if wall else 'NA'))
put(pd.DataFrame(resources),'resource_usage.tsv')
(R/'runtime_authorization.json').write_text(json.dumps({'user_instruction':'没有必要等，直接算就可以了','daytime_compute_authorized':True,'effective_scope':'this task; no waiting required','cpu_affinity_max':4,'BLAS_OpenMP_num_threads':1,'no_upstream_reanalysis':True},ensure_ascii=False,indent=2))
print('DELIVERY TABLES OK',flush=True)

a=pd.read_csv(R/'tables/association_primary.tsv',sep='\t')
a['evaluability']=np.where(a.P.notna(),'partially_evaluable_within_section','insufficient_for_inference')
a['association_conclusion']=np.where(a.P.isna(),'information_insufficient_for_test',np.where(a.BH_q<.05,'spatial_association_supported','no_clear_support_after_BH'))
put(a,'endpoint_evaluability.tsv')
