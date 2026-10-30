#!/usr/bin/env python3
"""Complete-section derived contrasts and neutral all-lineage displays; no inference."""
from summary_common import *
import argparse,time
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from circle_renderer import render_triplet
p=argparse.ArgumentParser();p.add_argument('method',choices=['COMMOT','SpatialCellChat']);p.add_argument('section');p.add_argument('--force',action='store_true');p.add_argument('--tables-only',action='store_true');a=p.parse_args()
method,section=a.method,a.section;out,source,variant=job_dirs(method,section);out.mkdir(parents=True,exist_ok=True);marker=out/'POSTPROCESS_COMPLETED.json'
if marker.exists() and not a.force and not a.tables_only:print('POSTPROCESS REUSED',method,section,flush=True);exit()
start=time.time();inp,info,meta,support,crosswalk=section_info(section);support=support.copy();support['type_evaluable']=support.type_evaluable.map(truth)
qf=read(inp/'Q_fine.tsv').set_index('spot_id');qm=read(inp/'Q_major.tsv').set_index('spot_id');coords=read(inp/'coords.tsv');assert list(qf.index)==meta.spot_id.tolist()==coords.spot_id.tolist()
labelmap=dict(zip(crosswalk.dropna(subset=['analysis_fine_label']).analysis_fine_label,crosswalk.dropna(subset=['analysis_fine_label']).major));labelmap.update({n:'Epithelial' for n in EPI})
hi,lo=EPI[:2];ss=support[(support.level=='fine')].set_index('node');state_ok={s:bool(ss.loc[s,'type_evaluable']) for s in [hi,lo]};hl_ok=all(state_ok.values());hl_status='evaluable' if hl_ok else 'network_only'
# An incomplete frozen score may be contrasted, but never enters complete20 primary summaries.
score_status='complete20' if int(info.get('n_mp10_genes_measured',0))==20 else 'score_coverage_incomplete'
metric='abundance_weighted_mean_cross_spot_strength' if method=='COMMOT' else 'total_mean_strength'
metric_label='G: abundance-weighted cross-spot mean' if method=='COMMOT' else 'Total mean strength: sum of LR F/K'
old=read(P/'mp10_subtype_screen_round3/tables/node_order.tsv');colors=dict(zip(old.node,old.color));colors.update({'Mural':'#a87c55',hi:'#d73027',lo:'#4575b4','T_NK_unresolved':'#999999'})
for r in crosswalk.itertuples():
 if pd.notna(r.analysis_fine_label) and r.analysis_fine_label not in colors:colors[r.analysis_fine_label]=colors.get(r.major,'#999999')
if method=='COMMOT':
 pinfo=json.loads((source/'parameters.json').read_text());lrn=int(pinfo['n_lr']);gene_n=int(pinfo['n_measured_genes'])
else:
 lrdata=read(source/'LR_selected.tsv');lrn=len(lrdata);gene_n=int(info.get('n_measured_genes',len(read(inp/'genes.tsv'))))
base=dict(method=method,section_id=section,dataset=unknown(info.get('dataset')),patient_uid=info.get('patient_uid',np.nan),patient_id_status=unknown(info.get('patient_id_status')),analysis_unit_uid=info.get('analysis_unit_uid',section),analysis_unit_type=unknown(info.get('analysis_unit_type')),tissue_condition=unknown(info.get('tissue_condition')),treatment_status=unknown(info.get('treatment_status')),score_coverage_status=score_status,input_score_coverage_status=info.get('score_coverage_status'),n_mp10_genes_measured=int(info.get('n_mp10_genes_measured',0)),hl_status=hl_status,n_spots=len(meta),n_epi_anchors=int(info.get('n_epi_anchors',meta.epi_anchor.map(truth).sum())),n_epi_high=int((meta.epi_state=='high').sum()),n_epi_low=int((meta.epi_state=='low').sum()),n_LR=lrn,n_measured_genes=gene_n,metric=metric,source_cache=str(source),variant=variant)
base.update({k:info.get(k,'unknown') for k in ['annotation_revision','tissue_evidence_level','tissue_evidence','tissue_evidence_file','clinical_annotation_source_sheet','clinical_annotation_source_rows','sample_disease_subtype','clinical_disease_background','sample_pathologic_stage','specimen_site','tumor_involvement_verified','annotation_limitations']})
rows=[];figures=[];checks=[]
for level,Q in [('major',qm),('fine',qf)]:
 e=read(source/f'{level}_edges.tsv').rename(columns={'total_mean_strength':'strength','evaluable':'computed','same_spot_flow_excluded':'same_spot_flow','summed_LR_supported_pair_denominators':'summed_LR_K'})
 nodes=list(Q.columns);targets=[n for n in nodes if n not in EPI];sx=support[support.level==level].set_index('node');eligible={n:bool(sx.loc[n,'type_evaluable']) for n in nodes}
 M=e.pivot(index='source',columns='target',values='strength').reindex(index=nodes,columns=nodes);uneligible=[n for n in nodes if not eligible[n]];M.loc[uneligible,:]=np.nan;M.loc[:,uneligible]=np.nan
 assert (M.stack().dropna()>=0).all();checks.append(dict(check=level+'_nonnegative_or_NE',pass_check=True))
 for direction in ['epi_to_TME','TME_to_epi']:
  low=M.loc[lo,targets] if direction=='epi_to_TME' else M.loc[targets,lo];high=M.loc[hi,targets] if direction=='epi_to_TME' else M.loc[targets,hi]
  paired=low.notna()&high.notna()&hl_ok;common=[n for n in targets if paired[n]];denlow=low[common].sum() if common else np.nan;denhigh=high[common].sum() if common else np.nan
  for n in targets:
   lin=n if level=='major' else labelmap[n];lnodes=[k for k in common if (k if level=='major' else labelmap[k])==lin];dl=low[lnodes].sum() if lnodes else np.nan;dh=high[lnodes].sum() if lnodes else np.nan
   ok=bool(paired[n]);vl=float(low[n]) if ok else np.nan;vh=float(high[n]) if ok else np.nan
   reason='NA' if ok else ('insufficient_Epi_high_or_low_support' if not hl_ok else 'target_type_support_insufficient' if not eligible[n] else 'aggregation_denominator_unavailable')
   zero='NE' if not ok else 'both_zero' if vl==0 and vh==0 else 'low_zero_only' if vl==0 else 'high_zero_only' if vh==0 else 'both_positive'
   row={**base,'network_level':level,'lineage':lin,'node':n,'computed':ok,'NE_reason':reason,'direction':direction,'strength_low':vl,'strength_high':vh,'delta_strength':vh-vl,'log2_HL':np.log2(vh/vl) if zero=='both_positive' else np.nan,'zero_category':zero,'effect_direction':'NE' if not ok else 'higher' if vh>vl else 'lower' if vh<vl else 'equal','share_all_TME_low':vl/denlow if denlow>0 else np.nan,'share_all_TME_high':vh/denhigh if denhigh>0 else np.nan,'share_within_lineage_low':vl/dl if dl>0 else np.nan,'share_within_lineage_high':vh/dh if dh>0 else np.nan,'all_TME_denominator_low':denlow,'all_TME_denominator_high':denhigh,'all_TME_denominator_labels':';'.join(common),'lineage_denominator_low':dl,'lineage_denominator_high':dh,'lineage_denominator_labels':';'.join(lnodes),'lineage_n_computed':len(lnodes),'target_support_spots':sx.loc[n,'n_positive_Q'],'target_abundance_support_spots':sx.loc[n,'n_abundance_ge_0_01'],'target_expected_abundance':sx.loc[n,'total_expected_abundance'],'target_effective_Q_mass':sx.loc[n,'effective_spot_mass'],'target_n_dominant_spots':sx.loc[n,'n_dominant_spots'],'identifiability_status':sx.loc[n,'identifiability_status']}
   for state,ename in [('low',lo),('high',hi)]:
    es,et=(ename,n) if direction=='epi_to_TME' else (n,ename);er=e[(e.source==es)&(e.target==et)].iloc[0]
    for col in ['total_flow','denominator','summed_LR_K','same_spot_flow','positive_cross_links','source_support_spots','target_support_spots']:
     if col in er:row[f'{col}_{state}']=er[col] if ok else np.nan
    row[f'epi_expected_abundance_{state}']=sx.loc[ename,'total_expected_abundance'];row[f'epi_effective_Q_mass_{state}']=sx.loc[ename,'effective_spot_mass'];row[f'epi_abundance_support_spots_{state}']=sx.loc[ename,'n_abundance_ge_0_01']
   row['delta_flow']=row.get('total_flow_high',np.nan)-row.get('total_flow_low',np.nan);row['mean_flow_direction_agrees']=(np.sign(row['delta_strength'])==np.sign(row['delta_flow'])) if ok else pd.NA
   rows.append(row)
 # All lineages treated identically; invalid full-section contrast does not produce empty circles.
 if hl_ok and not a.tables_only:
  groups=[('all_major',targets)] if level=='major' else [(lin,[n for n in targets if labelmap[n]==lin]) for lin in LINEAGES]
  for lineage,target_nodes in groups:
   if level=='fine' and lineage=='Endothelial':target_nodes=target_nodes+['Tip EC','Cycling EC']
   ns=[hi,lo]+target_nodes;matrix=M.reindex(index=ns,columns=ns).to_numpy();caption='All included protein LR, 3 classes; c2l proportion-weighted cross-spot attribution from mixed expression.\n'+('Incomplete MP10 panel: supplemental contrast; ' if score_status!='complete20' else '')+'NE = unavailable; scales differ across panels; no membrane-contact claim.'
   title=f'{method} | {section}\nVisium spot | '+('major' if level=='major' else 'fine: '+lineage)
   fr=render_triplet(matrix,ns,[colors.get(n,colors.get(lineage,'#999999')) for n in ns],hi,lo,out/'figures'/f'{level}_{safe(lineage)}',title,subtitle=caption,computed=[eligible.get(n,False) for n in ns],preview=section in ['GSE171351_A1','Inhousedata_T4469479','GSE171351_B1','GSE238145_B10'] and lineage in ['all_major','T/NK','Endothelial','Fibroblast','Mural'],metric=metric_label,epi_caption='Epi_H/L = epithelial components in MP10-high/low-associated anchor spots; shared TME Q')
   for f in fr:figures.append({**base,'network_level':level,'lineage':lineage,'figure_type':'circle','direction':'Epi to/from TME, both directions','source_data':str(source/f'{level}_edges.tsv'),'scale_definition':'same linear strength mapping across High/Low/Both; no positive edges hidden',**f})
C=pd.DataFrame(rows)
for sh in ['all_TME','within_lineage']:C['delta_share_'+sh]=C['share_'+sh+'_high']-C['share_'+sh+'_low']
for col in ['strength_low','strength_high','delta_strength','delta_share_within_lineage']:
 C['rank_'+col]=C.groupby(['network_level','lineage','direction'])[col].rank(method='min',ascending=False,na_option='keep')
save(C,out/'epi_state_type_contrasts.tsv');save(C[[c for c in C if 'share' in c or c in ['method','section_id','dataset','patient_uid','patient_id_status','tissue_condition','treatment_status','score_coverage_status','network_level','lineage','node','direction','computed','NE_reason','all_TME_denominator_labels','lineage_denominator_labels','lineage_n_computed']]],out/'lineage_shares.tsv')
# Unmodeled labels stay explicit NE and are never invented as zero-Q columns.
u=crosswalk[crosswalk.identifiability_status.eq('NE_NOT_IN_SPATIAL_REFERENCE')].copy();u['method']=method;u['computed']=False;u['strength_low']=np.nan;u['strength_high']=np.nan;save(u,out/'unavailable_types.tsv')
if not a.tables_only:
 if hl_ok:
  for level,lineage in [('major','all_major')]+[('fine',lin) for lin in LINEAGES]:
   g=C[(C.network_level==level)&((C.lineage==lineage) if level=='fine' else True)];names=list(g.node.drop_duplicates());fig,axs=plt.subplots(2,4,figsize=(17,max(7,.6*len(names)+3)),squeeze=False)
   for ii,direction in enumerate(['epi_to_TME','TME_to_epi']):
    z=g[g.direction==direction].set_index('node').reindex(names);y=np.arange(len(z));labels=[n+(' [NE]' if not bool(z.loc[n,'computed']) else '') for n in names]
    ax=axs[ii,0];ax.hlines(y,z.strength_low,z.strength_high,color='#aaa');ax.scatter(z.strength_low,y,color='#4575b4',s=28,label='Low');ax.scatter(z.strength_high,y,color='#d73027',s=28,label='High');ax.set_yticks(y,labels);ax.set_title(direction+' | mean strength');ax.legend(fontsize=8)
    ax=axs[ii,1];ax.barh(y,z.delta_strength,color=np.where(z.delta_strength>=0,'#d73027','#4575b4'));ax.axvline(0,color='black',lw=.6);ax.set_title('High - Low')
    for jj,key in [(2,'all_TME'),(3,'within_lineage')]:
     ax=axs[ii,jj];ax.hlines(y,z['share_'+key+'_low']*100,z['share_'+key+'_high']*100,color='#aaa');ax.scatter(z['share_'+key+'_low']*100,y,color='#4575b4');ax.scatter(z['share_'+key+'_high']*100,y,color='#d73027');ax.set_title(key.replace('_',' ')+' (%)')
    for ax in axs[ii]:ax.invert_yaxis();ax.grid(axis='x',alpha=.15);ax.spines[['right','top']].set_visible(False)
    for ax in axs[ii,1:]:ax.set_yticks(y,[])
   fig.suptitle(f'{method} | {section} | {level}: {lineage}\n{metric_label}',fontsize=14);fig.text(.5,.012,'Composition-weighted cross-spot attribution; matched H/L eligibility; no patient-level test. Single-subtype lineage share=1 is structural.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.04,1,.91]);stem=out/'figures'/('comparison_'+level+'_'+safe(lineage));fig.savefig(str(stem)+'.pdf',transparent=True);fig.savefig(str(stem)+'.png',transparent=True,dpi=300);plt.close(fig)
   figures.append({**base,'network_level':level,'lineage':lineage,'figure_type':'comparison','direction':'both Epi to/from TME','source_data':str(out/'epi_state_type_contrasts.tsv'),'scale_definition':'linear numeric axes; distinct mean/flow/share units','pdf':str(stem)+'.pdf','png':str(stem)+'.png','metric':metric})
 else:
  reason='No valid H/L comparison: each Epi state needs >=10 spots with mean abundance >=0.01 and total abundance >1.'
  text=f'{method} | {section}\nNetwork only; H/L comparison NE\n\n{reason}\n\nHigh spots: {base["n_epi_high"]}; Low spots: {base["n_epi_low"]}; anchors: {base["n_epi_anchors"]}\nThe full measured-position network is retained. No pseudo H/L circle is drawn.'
  fig,ax=plt.subplots(figsize=(11,5));ax.axis('off');ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,fontsize=11,wrap=True);stem=out/'figures/network_only_status';stem.parent.mkdir(exist_ok=True);fig.savefig(str(stem)+'.pdf',transparent=True);fig.savefig(str(stem)+'.png',transparent=True,dpi=300);plt.close(fig)
  (out/'HL_STATUS.md').write_text(f'# {section} / {method}\n\nnetwork_only；H/L 对照为 NE。\n\n{reason}\n完整位置网络与支持信息仍保留。\n')
  figures.append({**base,'network_level':'all','lineage':'all','figure_type':'network_only_status','direction':'NE','source_data':str(inp/'type_support.tsv'),'scale_definition':'not applicable; state explanation only','pdf':str(stem)+'.pdf','png':str(stem)+'.png','metric':metric})
# Lightweight geometry/Q/quality diagnostics, kept separate from communication inference.
xy=coords[['x_um','y_um']].to_numpy();tree=cKDTree(xy);geom={}
for radius in [105.,200.]:
 counts=tree.query_ball_point(xy,r=radius+1e-6,return_length=True)-1;geom[f'n_reachable_pairs_{int(radius)}um']=int(counts.sum());geom[f'fraction_reachable_pairs_{int(radius)}um']=counts.sum()/(len(xy)*(len(xy)-1));geom[f'median_neighbors_{int(radius)}um']=float(np.median(counts))
geom['nominal_tissue_area_proxy_um2']=len(xy)*100**2*np.sqrt(3)/2;geom['median_observed_nearest_neighbor_um']=float(np.median(tree.query(xy,k=2)[0][:,1]));geom['distance_scale_is_platform_not_observed_nn']=True
for state in ['high','low','other']:
 m=meta.epi_state.eq(state)
 for field in ['library_counts','n_genes_detected','MP10_frozen','epi_fraction_frozen']:
  geom[f'{state}_{field}_median']=float(meta.loc[m,field].median()) if m.any() else np.nan
save([{**base,**geom}],out/'batch_geometry_diagnostics.tsv')
for (level,lineage,direction),g in C.groupby(['network_level','lineage','direction']):
 for state in ['low','high']:
  s=g[f'share_within_lineage_{state}'].sum(min_count=1);assert pd.isna(s) or abs(s-1)<1e-9
  checks.append(dict(check=f'{level}_{lineage}_{direction}_{state}_lineage_share_sum',pass_check=True,value=s))
if not a.tables_only:
 F=pd.DataFrame(figures);F['relative_pdf']=F.pdf.map(lambda x:str(Path(x).relative_to(O)));F['relative_png']=F.png.map(lambda x:str(Path(x).relative_to(O)));save(F,out/'figure_index.tsv')
 circles=F[F.figure_type=='circle'];assert len(circles)==(27 if hl_ok else 0)
 for key,g in circles.groupby(['network_level','lineage']):assert len(g)==3 and g.edge_max.nunique()==1 and not g.has_self_loops.any() and g.hidden_positive_edges.eq(0).all()
 checks.append(dict(check='circle_count_and_triplet_scale',pass_check=True,value=len(circles)));save(checks,out/'postprocess_checks.tsv')
 atom_json({**base,'status':'completed','hl_status':hl_status,'figures':len(figures),'circles':len(circles),'n_contrasts':len(C),'n_evaluable_contrasts':int(C.computed.sum()),'elapsed_seconds':time.time()-start,'postprocess_contract':'all_sections_v2_1'},marker)
print('POSTPROCESS DONE',method,section,hl_status,'evaluable',int(C.computed.sum()),'figures',len(figures),flush=True)
