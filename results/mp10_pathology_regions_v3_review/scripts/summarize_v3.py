from run_v3 import *
import html,platform,re
T=O/'tables';F=O/'figures';T.mkdir(exist_ok=True);F.mkdir(exist_ok=True)
identity=pd.read_csv(R/'spatial_ecology_top20_v2/patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv',sep='\t').set_index('section_id')
summary=[];allspots=[];domains=[];rois=[];bands=[];changes=[];reviews=[];biomarkers=[];contactsumm=[];localreview=[];cohortcontacts=[]
features={'MP10':'MP10_original','malignant_fraction':'malignant_fraction','Epithelial':'c2l_fraction_Epithelial','Treg':'c2l_fraction_Treg','EC':'c2l_Endothelial','CAF':'c2l_stroma'}
def composition(z):
 z=z.copy();z['c2l_Endothelial']=z[[f'c2l_fraction_{s} EC' for s in ['Arterial','Capillary','Venous','Lymphatic']]].sum(axis=1);z['c2l_stroma']=z[[f'c2l_fraction_{s}' for s in ['Matrix_CAF','iCAF','rCAF','vCAF']]].sum(axis=1);return z
for sec,m in M.iterrows():
 p=O/'by_dataset'/m.dataset_id/sec
 if not (p/'complete.json').exists():raise RuntimeError('Section not complete: '+sec)
 d=composition(pd.read_parquet(p/'spot_regions_v3.parquet'));assert len(d)==m.n_spots
 row=json.loads((p/'complete.json').read_text());co=json.loads((p/'coordinate_scale.json').read_text());row['distance_estimable_fraction']=d.signed_distance_um.notna().mean();row['roi_spot_footprint_um2']=row['n_inside']*10000*np.sqrt(3)/2;row['affine_residual_max_pitch']=co['affine_residual_max_pitch'];row['image_available']=co['image_available'];row['v2_tumor_outside_v3']=int((d.v2_tumor_mask&d.lesion_location.eq('outside')).sum());row['v3_inside_non_v2_tumor']=int((~d.v2_tumor_mask&d.lesion_location.eq('inside')).sum());row['v3_inside_below_old_epi_gate']=int((d.lesion_location.eq('inside')&(d.epi_fraction_gate<.25)).sum());row['v2_gray_total']=int(d.v2_tissue_region.isin(['Other tissue','Uncertain/conflicting']).sum());row['boundary_both_sides_20spots']=bool(d.distance_band.eq('Inner margin').sum()>=20 and d.distance_band.eq('Outer margin').sum()>=20);row['interior_outside_20spots']=bool(d.distance_band.eq('Tumor interior').sum()>=20 and d.distance_band.eq('Outside').sum()>=20);row['interpretation_caution']=sec in ['GSE285715_N2','GSE285715_N3','GSE285715_T1','GSE285715_T2','GSE285715_T3','CanDis_2026_Sample8'] or row['image_alignment_caution']
 ident=identity.loc[sec];meta={k:ident[k] for k in ['patient_uid','patient_id_status','analysis_unit_uid','analysis_unit_type']};row.update(meta);summary.append(row)
 for k,v in meta.items():d[k]=v
 allspots.append(d);dt=pd.read_csv(p/'domain_annotations.tsv',sep='\t');domains.append(dt)
 ed=pd.read_csv(p/'lesion_boundary_edges.tsv.gz',sep='\t')
 for roi,z in d[d.lesion_location.eq('inside')].groupby('lesion_roi_id'):
  ee=ed[ed.lesion_roi_id==roi];flags=[]
  if len(z)<10:flags.append('small candidate retained; inspect local contour')
  if z.signed_distance_um.notna().mean()<.25:flags.append('depth mostly censored; use location only outside observed edge neighborhood')
  if sec in ['GSE285715_T1','GSE285715_T2','GSE285715_T3']:flags.append('weak/low-confidence CNA or expression-signature dependence')
  rois.append({'section_id':sec,'dataset_id':m.dataset_id,'lesion_roi_id':roi,'n_spots':len(z),'n_v2_non_tumor_inside':int((~z.v2_tumor_mask).sum()),'n_valid_outer_edges':int(ee.boundary_type.eq('biological_lesion_outer').sum()),'distance_estimable_fraction':z.signed_distance_um.notna().mean(),'malignant_fraction_median':z.malignant_fraction.median(),'review_status':C['roi_review_status'],'review_flags':'; '.join(flags)})
  if flags:localreview.append({'section_id':sec,'local_object':roi,'x_um_min':z.x_um.min(),'x_um_max':z.x_um.max(),'y_um_min':z.y_um.min(),'y_um_max':z.y_um.max(),'issue':'; '.join(flags)})
 for w in C['sensitivity_half_widths_um']:
  for label,n in d[f'distance_band_{w}um'].value_counts().items():bands.append({'section_id':sec,'dataset_id':m.dataset_id,'half_width_um':w,'distance_band':label,'n_spots':n})
 for (old,loc,band),z in d.groupby(['v2_tissue_region','lesion_location','main_region'],dropna=False):changes.append({'section_id':sec,'dataset_id':m.dataset_id,'v2_identity':old,'v3_location':loc,'v3_main_region':band,'n_spots':len(z)})
 for band,z in d.groupby('distance_band'):
  a={'section_id':sec,'dataset_id':m.dataset_id,**meta,'boundary_scope':'lesion_outer','distance_band':band,'n_spots':len(z),'n_distance_estimable':z.signed_distance_um.notna().sum()}
  for key,col in features.items():a[key+'_median']=z[col].median();a[key+'_mean']=z[col].mean();a[key+'_q25']=z[col].quantile(.25);a[key+'_q75']=z[col].quantile(.75)
  biomarkers.append(a)
 # Tumor-nest results remain separate, restrict follow-up summary to contacts within candidate lesion ROI.
 ct=pd.read_csv(p/'tumor_nest_contact_spots.tsv.gz',sep='\t');assert ct.barcode.tolist()==d.barcode.tolist();cohortcontacts.append(ct);zz=d.copy();zz['nest_distance']=ct.signed_distance_um.to_numpy();zz['contact_status']=ct.distance_status.to_numpy();zz=zz[zz.lesion_location.eq('inside')&zz.contact_status.eq('observed_contact_distance')]
 zz['contact_band']=pd.cut(zz.nest_distance,[-np.inf,-200,0,200,np.inf],labels=['Nest interior','Nest inner contact','TME outer contact','TME distant'],right=False)
 for label,z in zz.groupby('contact_band',observed=True):
  a={'section_id':sec,'dataset_id':m.dataset_id,**meta,'boundary_scope':'tumor_nest_contact','contact_band':str(label),'n_spots':len(z)}
  for key,col in features.items():a[key+'_median']=z[col].median()
  contactsumm.append(a)
 if row['n_location_unresolved']:
  z=d[d.lesion_location=='unlocalized'];localreview.append({'section_id':sec,'local_object':'unlocalized spots','x_um_min':z.x_um.min(),'x_um_max':z.x_um.max(),'y_um_min':z.y_um.min(),'y_um_max':z.y_um.max(),'issue':'source/morphology/malignancy conflict; inspect explicit gray spots, not entire section'})
 if row['image_alignment_caution']:localreview.append({'section_id':sec,'local_object':'image transform','issue':f"max affine residual = {co['affine_residual_max_pitch']:.3f} pitch; nominal array geometry retained, H&E overlay requires local review"})
 # Compact standalone main figure for direct review (all sections, including no-lesion controls).
 fig,ax=plt.subplots(figsize=(6,6));ax.scatter(d.x,d.y,c=[PAL[x] for x in d.main_region],s=6,marker='h',linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(sec+'\nLesion space | w = 200 um',fontsize=10);ax.legend(handles=[Patch(color=v,label=k) for k,v in PAL.items() if k in set(d.main_region)],loc='upper center',bbox_to_anchor=(.5,0),ncol=2,frameon=False,fontsize=7);fig.subplots_adjust(bottom=.19);fig.savefig(p/'05_main_regions.png',dpi=300);fig.savefig(p/'05_main_regions.pdf');plt.close(fig)
S=pd.DataFrame(summary);D=pd.concat(allspots,ignore_index=True);assert not D.duplicated(['dataset_id','section_id','barcode']).any();S.to_csv(T/'section_summary.tsv',sep='\t',index=False);pd.concat(domains,ignore_index=True).to_csv(T/'domain_annotations.tsv.gz',sep='\t',index=False);pd.DataFrame(rois).to_csv(T/'lesion_roi_summary.tsv',sep='\t',index=False);pd.DataFrame(localreview).to_csv(T/'local_review_targets.tsv',sep='\t',index=False);pd.DataFrame(bands).to_csv(T/'bandwidth_sensitivity.tsv',sep='\t',index=False);pd.DataFrame(changes).to_csv(T/'v2_to_v3_transitions.tsv',sep='\t',index=False);D.to_parquet(T/'all_spots_v3.parquet',index=False)
core=['dataset_id','section_id','barcode','x','y','x_um','y_um','spatial_domain','domain_component_id','domain_annotation','lesion_roi_id','lesion_location','boundary_scope','reference_lesion_roi_id','signed_distance_um','distance_status','distance_band','main_region','location_confidence','malignant_fraction','SpaCET_signature_source','CNV_support','malignant_evidence_confidence','quality_status','assignment_basis','roi_review_status','MP10_original','patient_uid','analysis_unit_uid']
D[core].to_csv(T/'all_spot_labels_v3.tsv.gz',sep='\t',index=False);pd.concat(cohortcontacts,ignore_index=True).to_parquet(T/'tumor_nest_contact_spots.parquet',index=False)
B=pd.DataFrame(biomarkers);B.to_csv(T/'lesion_outer_region_MP10_TME_summary.tsv',sep='\t',index=False);pd.DataFrame(contactsumm).to_csv(T/'tumor_nest_contact_MP10_TME_summary.tsv',sep='\t',index=False)
# Section paired contrasts: no spot-wise significance tests; repeated patient sections averaged at the verified analysis unit.
contrasts=[]
for sec,z in B.groupby('section_id'):
 z=z.set_index('distance_band')
 for a,b,name in [('Tumor interior','Outside','interior_minus_outside'),('Inner margin','Outer margin','inner_minus_outer_margin')]:
  if a not in z.index or b not in z.index or min(z.loc[a,'n_spots'],z.loc[b,'n_spots'])<20:continue
  meta=z.loc[a];q={'section_id':sec,'dataset_id':meta.dataset_id,'patient_uid':meta.patient_uid,'analysis_unit_uid':meta.analysis_unit_uid,'analysis_unit_type':meta.analysis_unit_type,'contrast':name,'n_a':int(meta.n_spots),'n_b':int(z.loc[b,'n_spots'])}
  for key in features:q[key+'_median_difference']=meta[key+'_median']-z.loc[b,key+'_median']
  contrasts.append(q)
K=pd.DataFrame(contrasts);K.to_csv(T/'section_paired_contrasts.tsv',sep='\t',index=False)
cols=[x for x in K if x.endswith('_difference')];U=K.groupby(['dataset_id','analysis_unit_uid','analysis_unit_type','contrast'],dropna=False)[cols].mean().reset_index();U.to_csv(T/'analysis_unit_paired_contrasts.tsv',sep='\t',index=False)
cons=[]
for (ds,cc),z in U.groupby(['dataset_id','contrast']):
 for key in features:
  v=z[key+'_median_difference'];cons.append({'dataset_id':ds,'contrast':cc,'feature':key,'n_analysis_units':len(v),'n_positive':int((v>0).sum()),'n_negative':int((v<0).sum()),'median_difference':v.median(),'interpretation':'descriptive direction; no spot-level independent-replicate inference'})
pd.DataFrame(cons).to_csv(T/'dataset_direction_consistency.tsv',sep='\t',index=False)
# Purity-stratified descriptive check uses original scores, not a fitted cell-state correction.
purity=[]
for sec,z in D[D.distance_status.eq('estimable')].groupby('section_id'):
 z=z.copy();z['purity_bin']=pd.cut(z.malignant_fraction,[-.001,.25,.5,.75,1.001],labels=['0-.25','.25-.5','.5-.75','.75-1'])
 for (pb,band),zz in z.groupby(['purity_bin','distance_band'],observed=True):purity.append({'section_id':sec,'dataset_id':zz.dataset_id.iloc[0],'malignant_fraction_bin':str(pb),'distance_band':band,'n_spots':len(zz),'MP10_median':zz.MP10_original.median(),'note':'shared model-estimated fraction stratum; not purified tumor-cell state'})
pd.DataFrame(purity).to_csv(T/'MP10_malignant_fraction_stratified.tsv',sep='\t',index=False)
# Per-dataset all-section atlas.
for ds,ss in S.groupby('dataset_id'):
 n=len(ss);nc=min(4,n);nr=int(np.ceil(n/nc));fig,axs=plt.subplots(nr,nc,figsize=(nc*3.6,nr*3.5+1),squeeze=False)
 for ax,sec in zip(axs.ravel(),ss.section_id):
  z=D[D.section_id==sec];ax.scatter(z.x,z.y,c=[PAL[x] for x in z.main_region],s=2,marker='h',linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(sec.replace(ds+'_', ''),fontsize=8)
 for ax in axs.ravel()[n:]:ax.axis('off')
 fig.legend(handles=[Patch(color=v,label=k) for k,v in PAL.items()],loc='lower center',ncol=4,fontsize=7,frameon=False);fig.suptitle(ds+' | all sections | fixed 200 um half-width',fontsize=12);fig.subplots_adjust(bottom=.13,top=.92,hspace=.3);fig.savefig(F/(ds+'_overview.png'),dpi=300);fig.savefig(F/(ds+'_overview.pdf'));plt.close(fig)
fig,axs=plt.subplots(1,4,figsize=(14,4));dsnames=sorted(S.dataset_id.unique());dc={x:plt.get_cmap('tab10')(i) for i,x in enumerate(dsnames)}
for ax,key in zip(axs,['MP10','Treg','EC','CAF']):
 z=K[K.contrast=='interior_minus_outside'];
 for i,ds in enumerate(dsnames):
  v=z[z.dataset_id==ds][key+'_median_difference'].to_numpy();jitter=np.linspace(-.15,.15,len(v)) if len(v)>1 else np.zeros(len(v));ax.scatter(i+jitter,v,s=18,color=dc[ds],alpha=.8)
 ax.axhline(0,c='gray',lw=.7);ax.set_xticks(range(len(dsnames)),[x.replace('__UNVERIFIED','') for x in dsnames],rotation=80,fontsize=7);ax.set_title(key);ax.set_ylabel('Interior − outside; section median difference',fontsize=8)
fig.suptitle('Descriptive section contrasts; ≥20 estimable spots per region\nComposition contributed to ROI annotation: not independent validation',fontsize=10);fig.tight_layout();fig.savefig(F/'MP10_TME_section_contrasts.png',dpi=300);fig.savefig(F/'MP10_TME_section_contrasts.pdf');plt.close(fig)
# Actual Cottrazm diagnostic outputs, same coordinates; agreement is not accuracy.
comp=[]
for sec in C['cottrazm_sections']:
 cp=O/'comparators'/sec;lab=cp/'cottrazm_labels.tsv';z=D[D.section_id==sec]
 if not lab.exists():continue
 a=pd.read_csv(lab,sep='\t').set_index('barcode').reindex(z.barcode);assert a.Cottrazm_region.notna().all();pal={'Malignant':'#bd3345','Boundary':'#f5a14d','Non-malignant':'#4389b6','Unassigned':'#999999'}
 for (v3,ct),nn in pd.DataFrame({'v3':z.main_region.to_numpy(),'ct':a.Cottrazm_region.to_numpy()}).value_counts().items():comp.append({'section_id':sec,'v3_main_region':v3,'Cottrazm_adapted_module_region':ct,'n_spots':nn})
 fig,axs=plt.subplots(1,3,figsize=(15,5));vals=[z.v2_tissue_region,z.main_region,a.Cottrazm_region];pals=[{'Tumor':'#bd3345','Normal urothelium':'#e2b637','Non-tumor stromal region':'#479b8d','Immune-rich region':'#477dc1','Smooth-muscle-rich region':'#8b67ad','Other tissue':'#d3c6b6','Uncertain/conflicting':'#999999'},PAL,pal]
 for ax,val,pa,title in zip(axs,vals,pals,['v2 identity','v3 lesion space','Cottrazm BoundaryDefine (adapted inputs)']):
  ax.scatter(z.x,z.y,c=[pa.get(x,'#999999') for x in val],s=3,marker='h',linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off');ax.set_title(title,fontsize=10);ax.legend(handles=[Patch(color=v,label=k) for k,v in pa.items() if k in set(val)],loc='upper center',bbox_to_anchor=(.5,0),fontsize=6,frameon=False,ncol=2)
 fig.suptitle(sec+' | diagnostic comparison; no independent pathology accuracy estimate');fig.subplots_adjust(bottom=.22,wspace=.25);fig.savefig(cp/'comparison.png',dpi=300);fig.savefig(cp/'comparison.pdf');plt.close(fig)
pd.DataFrame(comp).to_csv(T/'Cottrazm_diagnostic_crosswalk.tsv',sep='\t',index=False)
# Browser-friendly index, relative links make the whole result directory portable.
rows=[]
for r in S.itertuples():
 base=f'by_dataset/{r.dataset_id}/{r.section_id}';rows.append(f'<tr><td>{r.dataset_id}</td><td>{r.section_id}</td><td>{r.n_spots}</td><td>{r.n_rois}</td><td>{r.distance_estimable_fraction:.1%}</td><td><a href="{base}/01_v3_core.png">H&E / ROI / distance</a> · <a href="{base}/02_v2_SpaCET_v3_comparison.png">v2 / SpaCET / v3</a> · <a href="{base}/05_main_regions.png">main</a> · <a href="{base}/spot_regions_v3.tsv.gz">labels</a></td></tr>')
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>BLCA lesion regions v3</title><style>body{font:15px sans-serif;margin:30px}td,th{padding:8px;border-bottom:1px solid #ddd}a{color:#176da2}table{border-collapse:collapse}</style><h1>BLCA lesion regions v3</h1><p>Computational candidate ROI; no documented pathologist confirmation. Main outer boundary and secondary nest contacts are separate. Pink/pale blue preserve inside/outside where distance is unestimated.</p><p><a href="RESULTS.md">Results</a> · <a href="METHODS.md">Methods</a> · <a href="tables/section_summary.tsv">Section summary</a> · <a href="tables/all_spot_labels_v3.tsv.gz">All labels</a></p><table><tr><th>Dataset</th><th>Section</th><th>Spots</th><th>ROI</th><th>Distance estimable</th><th>Results</th></tr>'+''.join(rows)+'</table>')
S[['dataset_id','section_id','n_spots','n_rois','distance_estimable_fraction']].assign(result_path=[str(O/'by_dataset'/r.dataset_id/r.section_id) for r in S.itertuples()]).to_csv(T/'result_index.tsv',sep='\t',index=False)
versions={'Python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__,'shapely':__import__('shapely').__version__,'matplotlib':matplotlib.__version__,'seed':C['seed']};(O/'logs/python_versions.json').write_text(json.dumps(versions,indent=2))
print('COMPLETE',len(S),len(D),'ROIs',S.n_rois.sum(),'estimable',S.n_distance_estimable.sum(),flush=True)
print(S[['n_inside','n_location_unresolved','v2_gray_total','v2_gray_resolved_inside','v2_gray_resolved_outside','v3_inside_non_v2_tumor','v2_tumor_outside_v3']].sum().to_string(),flush=True)
