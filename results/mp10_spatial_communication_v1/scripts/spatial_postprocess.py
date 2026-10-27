from runtime import *
import numpy as np,pandas as pd,json,argparse,re
import matplotlib.pyplot as plt
from circle_renderer import render_triplet
p=argparse.ArgumentParser();p.add_argument('method',choices=['COMMOT','SpatialCellChat']);p.add_argument('section');p.add_argument('--variant',default=None);args=p.parse_args()
method=args.method;section=args.section;variant=args.variant or ('primary2dNN' if method=='COMMOT' else 'main');folder='commot' if method=='COMMOT' else 'spatialcellchat';out=O/folder/section/variant
metric='G: abundance-weighted cross-spot mean' if method=='COMMOT' else 'Total mean strength: sum of LR F/K'
hi='Epi_MP10_high';lo='Epi_MP10_low';ep=[hi,lo,'Epi_MP10_other']
read=lambda p:pd.read_csv(p,sep='\t')
cross=read(O/'input_review/spatial_label_crosswalk.tsv').query('section_id==@section')
labelmap=dict(zip(cross.dropna(subset=['analysis_fine_label']).analysis_fine_label,cross.dropna(subset=['analysis_fine_label']).major));labelmap.update({e:'Epithelial' for e in ep})
lineages=['B','Plasma','T/NK','Myeloid','Fibroblast','Endothelial','Mast','Mural']
old=read(P/'mp10_subtype_screen_round3/tables/node_order.tsv');colors=dict(zip(old.node,old.color));majorcolors={'Mural':'#a87c55','Epi_MP10_high':'#d73027','Epi_MP10_low':'#4575b4','T_NK_unresolved':'#999999'};colors.update(majorcolors)
for r in cross.itertuples():
 if pd.notna(r.analysis_fine_label) and r.analysis_fine_label not in colors:colors[r.analysis_fine_label]=colors.get(r.major,'#999999')
type_support=read(O/'input_review/type_support.tsv').query('section_id==@section');pilot=read(O/'input_review/pilot_sections.tsv').query('section_id==@section').iloc[0]
allrows=[];figrows=[];edgerows=[]
for level in ['major','fine']:
 E=read(out/f'{level}_edges.tsv');E=E.rename(columns={'total_mean_strength':'strength','evaluable':'computed','same_spot_flow_excluded':'same_spot_flow','summed_LR_supported_pair_denominators':'summed_LR_K'})
 M=E.pivot(index='source',columns='target',values='strength');nodes=list(read(O/f'input_review/sections/{section}/Q_{level}.tsv').columns[1:]);M=M.reindex(index=nodes,columns=nodes)
 targets=[n for n in nodes if n not in ep];sidx=type_support.query('level==@level').set_index('node')
 # Actual paired eligibility and zero denominator retained explicitly.
 eligible={n:bool(sidx.loc[n,'type_evaluable']) for n in nodes};invalid=[n for n in nodes if not eligible[n]];M.loc[invalid,:]=np.nan;M.loc[:,invalid]=np.nan
 for direction in ['epi_to_TME','TME_to_epi']:
  low=M.loc[lo,targets] if direction=='epi_to_TME' else M.loc[targets,lo]
  high=M.loc[hi,targets] if direction=='epi_to_TME' else M.loc[targets,hi]
  common=low.notna()&high.notna();denlow=low[common].sum();denhigh=high[common].sum()
  for n in targets:
   lin=n if level=='major' else labelmap[n];group=[g for g in targets if (g if level=='major' else labelmap[g])==lin and common[g]]
   dl=low[group].sum() if group else np.nan;dh=high[group].sum() if group else np.nan
   row=dict(method=method,section_id=section,dataset=pilot.dataset,patient_uid=pilot.patient_uid,patient_id_status=pilot.patient_id_status,analysis_unit_uid=pilot.analysis_unit_uid,variant=variant,network_level=level,lineage=lin,node=n,direction=direction,computed=bool(common[n]),strength_low=low[n],strength_high=high[n],delta_strength=high[n]-low[n],share_all_TME_low=low[n]/denlow if denlow>0 else np.nan,share_all_TME_high=high[n]/denhigh if denhigh>0 else np.nan,share_within_lineage_low=low[n]/dl if dl>0 else np.nan,share_within_lineage_high=high[n]/dh if dh>0 else np.nan,all_TME_denominator_low=denlow,all_TME_denominator_high=denhigh,lineage_denominator_low=dl,lineage_denominator_high=dh,lineage_labels=';'.join(group),lineage_n_computed=len(group),n_spots=int(pilot.n_spots),n_epi_high=int(pilot.n_epi_high),n_epi_low=int(pilot.n_epi_low),target_support_spots=sidx.loc[n,'n_positive_Q'],target_effective_spot_mass=sidx.loc[n,'effective_spot_mass'],target_expected_abundance=sidx.loc[n,'total_expected_abundance'],target_n_dominant_spots=sidx.loc[n,'n_dominant_spots'],identifiability_status=sidx.loc[n,'identifiability_status'],metric=metric)
   for state,en in [('low',lo),('high',hi)]:
    s,t=(en,n) if direction=='epi_to_TME' else (n,en);r=E[(E.source==s)&(E.target==t)].iloc[0]
    for col in ['total_flow','denominator','summed_LR_K','same_spot_flow','positive_cross_links']:
     if col in r:row[f'{col}_{state}']=r[col]
   allrows.append(row)
 E.insert(0,'section_id',section) if 'section_id' not in E else None;E['method']=method;E['network_level']=level;E['variant']=variant;edgerows.append(E)
 # Only primary gets main triplets; wide sensitivity tables and panels suffice.
 if variant in ['primary2dNN','main']:
  groups=[('major',targets)] if level=='major' else [(lin,[n for n in targets if labelmap[n]==lin]) for lin in ['T/NK','Endothelial','Fibroblast','Mural','B','Plasma','Myeloid','Mast']]
  for lin,gn in groups:
   # Explicitly acknowledge important EC types unavailable in accepted spatial reference.
   if level=='fine' and lin=='Endothelial':gn=gn+['Tip EC','Cycling EC']
   ns=[hi,lo]+gn;mm=M.reindex(index=ns,columns=ns).to_numpy();safe=re.sub(r'[^A-Za-z0-9_-]+','_',lin)
   subtitle='All included protein LR (3 classes); c2l proportion-weighted cross-spot attribution from mixed measured expression.\nTypes expressible in the reference only; scales differ across panels; neighbour contact does not prove membrane contact.'
   fs=render_triplet(mm,ns,[colors.get(n,colors.get(lin,'#999999')) for n in ns],hi,lo,out/'figures'/f'{level}_{safe}',f'{method} | {section}\nVisium spot | {level} {lin}',subtitle=subtitle,computed=[eligible.get(n,False) for n in ns],preview=lin in ['major','T/NK','Endothelial','Fibroblast'],metric=metric,epi_caption='Epi_H/L = epithelial components in MP10-high/low-associated anchor spots; shared TME composition')
   for f in fs:figrows.append(dict(method=method,section_id=section,variant=variant,network_level=level,lineage=lin,figure_type='circle',direction='both directions Epi to/from TME only',source_matrix=str(out/f'{level}_edges.tsv'),**f))
C=pd.DataFrame(allrows)
for share in ['all_TME','within_lineage']:C['delta_share_'+share]=C['share_'+share+'_high']-C['share_'+share+'_low']
for field in ['strength_low','strength_high','delta_strength','delta_share_within_lineage']:
 C['rank_'+field]=C.groupby(['network_level','lineage','direction'])[field].rank(method='min',ascending=False,na_option='keep')
C.to_csv(out/'epi_state_type_contrasts.tsv',sep='\t',index=False,na_rep='NA');pd.concat(edgerows).to_csv(out/'all_group_edges.tsv',sep='\t',index=False,na_rep='NA')
C[[c for c in C if 'share' in c or c in ['method','section_id','variant','network_level','lineage','node','direction','lineage_labels','lineage_n_computed','computed']]].to_csv(out/'lineage_shares.tsv',sep='\t',index=False,na_rep='NA')
# Focus panels, both directions and denominators visible, no patients/spots-as-replicates tests.
for lin in ['T/NK','Endothelial','Fibroblast','Mural']:
 sub=C[(C.network_level=='fine')&(C.lineage==lin)];names=list(sub.node.drop_duplicates());fig,axs=plt.subplots(2,4,figsize=(17,max(7,.6*len(names)+3)),squeeze=False)
 for row,d in enumerate(['epi_to_TME','TME_to_epi']):
  g=sub[sub.direction==d].set_index('node').reindex(names);y=np.arange(len(g));a=axs[row,0]
  a.hlines(y,g.strength_low,g.strength_high,color='#aaaaaa');a.scatter(g.strength_low,y,color='#4575b4',label='Low',s=30);a.scatter(g.strength_high,y,color='#d73027',label='High',s=30);a.set_title(d+' | strength');a.set_yticks(y,names);a.legend(fontsize=8)
  a=axs[row,1];a.barh(y,g.delta_strength,color=np.where(g.delta_strength>=0,'#d73027','#4575b4'));a.axvline(0,color='black',lw=.7);a.set_title('High - Low strength')
  for j,key in [(2,'all_TME'),(3,'within_lineage')]:
   a=axs[row,j];a.hlines(y,g['share_'+key+'_low']*100,g['share_'+key+'_high']*100,color='#aaa');a.scatter(g['share_'+key+'_low']*100,y,color='#4575b4');a.scatter(g['share_'+key+'_high']*100,y,color='#d73027');a.set_title(key.replace('_',' ')+' share (%)')
  for a in axs[row]:a.invert_yaxis();a.grid(axis='x',alpha=.15);a.spines[['right','top']].set_visible(False)
  for a in axs[row,1:]:a.set_yticks(y,[])
 fig.suptitle(f'{method} | {section} | {lin} | {variant}\n{metric}',fontsize=14);fig.text(.5,.01,'Proportion-weighted cross-spot attribution; mixed measured expression. Exploratory single-section comparison, no patient-level test.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.04,1,.91]);stem=out/'figures'/('comparison_'+re.sub(r'[^A-Za-z0-9_-]+','_',lin));stem.parent.mkdir(exist_ok=True);fig.savefig(str(stem)+'.pdf',transparent=True);fig.savefig(str(stem)+'.png',transparent=True,dpi=300)
 if lin=='T/NK':fig.savefig(str(stem)+'_white.png',facecolor='white',dpi=130)
 plt.close(fig);figrows.append(dict(method=method,section_id=section,variant=variant,network_level='fine',lineage=lin,figure_type='comparison',pdf=str(stem)+'.pdf',png=str(stem)+'.png',metric=metric))
if variant in ['primary2dNN','main']:
 inp=O/'input_review/sections'/section;meta=read(inp/'meta.tsv');xy=read(inp/'coords.tsv');q=read(inp/'Q_fine.tsv');sr=read(out/('spot_communication.tsv' if method=='COMMOT' else 'spot_sender_receiver.tsv'))
 keys=['MP10_frozen','Treg','Capillary EC','Matrix_CAF','cross_sender','cross_receiver'];values=[meta.MP10_frozen,q.Treg,q['Capillary EC'],q.Matrix_CAF,sr['cross_sender' if method=='COMMOT' else 'outgoing_cross_spot'],sr['cross_receiver' if method=='COMMOT' else 'incoming_cross_spot']]
 fig,axs=plt.subplots(2,3,figsize=(14,9))
 for ax,k,v in zip(axs.flat,keys,values):
  sc=ax.scatter(xy.x_um,xy.y_um,c=v,s=5,cmap='viridis',rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(k);ax.set_xlabel('array x (um)');ax.set_ylabel('array y (um)');fig.colorbar(sc,ax=ax,shrink=.7)
 fig.suptitle(f'{method} | {section} | complete Visium section');fig.tight_layout();stem=out/'figures/spatial_supplement';fig.savefig(str(stem)+'.pdf',transparent=True);fig.savefig(str(stem)+'.png',transparent=True,dpi=300);plt.close(fig);figrows.append(dict(method=method,section_id=section,variant=variant,network_level='spot',lineage='all',figure_type='spatial',pdf=str(stem)+'.pdf',png=str(stem)+'.png',metric='continuous MP10 / Q / total cross-spot flow'))
pd.DataFrame(figrows).to_csv(out/'figure_index.tsv',sep='\t',index=False);print('POSTPROCESS DONE',method,section,variant,'figures',len(figrows),flush=True)
