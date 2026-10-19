from core import *
from stats import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(fig,p):
 guard('review:'+p.name,8*1024**2);fig.savefig(p.with_suffix('.png'),dpi=300,transparent=True,bbox_inches='tight');fig.savefig(p.with_suffix('.pdf'),transparent=True,bbox_inches='tight');plt.close(fig)

def choose(p):
 hits=p[p.q_global_BH.lt(.05)];level='GLOBAL_BH_HIT'
 if not len(hits):hits=p[p.q_within_BH.lt(.05)];level='WITHIN_COMBO_EXPLORATORY'
 if not len(hits):hits=p[p.p.notna()];level='NON_SIGNIFICANT_INFORMATION'
 if not len(hits):return hits,level
 hits=hits.copy();hits['inhouse_priority']=hits.contains_inhouse & hits.inhouse_actual_support
 summ=pd.read_csv(OUT/'combination_summary.tsv',sep='\t',dtype={'combo_id':str});hits=hits.merge(summ[['combo_id','stability_min']],on='combo_id');hits=hits.sort_values(['inhouse_priority','n_units','stability_min','q_global_BH','q_within_BH'],ascending=[False,False,False,True,True]);selected=[];combos=[];seen=set()
 # Composition-dominant label avoids three nearly identical examples when available.
 for r in hits.itertuples():
  if r.combo_id in combos:continue
  prof=pd.read_parquet(OUT/'combinations'/r.combo_id/'profiles.parquet');k=int(r.target_id.replace('niche',''));v=prof[prof.niche.eq(k)].groupby('celltype')['mean'].mean();signature=tuple(v.nlargest(3).index)
  if signature in seen and len(set(hits.combo_id)-set(combos))>3:continue
  combos.append(r.combo_id);seen.add(signature)
  if len(combos)==3:break
 for cid in hits.combo_id:
  if len(combos)==3:break
  if cid not in combos:combos.append(cid)
 for cid in combos:selected.append(hits[hits.combo_id.eq(cid)].head(3))
 return pd.concat(selected,ignore_index=True),level

def run():
 c=pd.read_csv(OUT/'combination_summary.tsv',sep='\t',dtype={'combo_id':str});assert len(c)==99 and c.status.eq('COMPLETE').all(),'REVIEWS_REQUIRE_COMPLETE_MAIN_SCAN'
 p=pd.read_parquet(OUT/'primary_relations.parquet');assert p.global_status.eq('FINAL').all();selected,level=choose(p);table(OUT/'review/selected_relations.tsv',selected);js(OUT/'review/selection.json',dict(level=level,rule='global hits then within hits then informative non-significant; inhouse estimable+same direction, N, structure, q, composition diversity; <=3 combos x <=3 relations',combos=selected.combo_id.unique().tolist() if len(selected) else []))
 fig,ax=plt.subplots(figsize=(11,4));xx=np.arange(99);ax.scatter(xx,c.global_hits.fillna(0),c=np.where(c.contains_inhouse,'#cc5500','#3f6ba4'),s=18);ax.set(xlabel='Fixed execution order (all seven first)',ylabel='Global BH primary relations',title='99-combination exhaustive search; shared data, not independent validation');save(fig,OUT/'review/global_overview')
 from spotted import section,source,effects
 import importlib,rings
 ringrun=importlib.reload(rings).run
 master=pd.read_parquet(OUT/'shared/master_anchors.parquet');base=np.load(OUT/'shared/base_domains.npz')['base_score'];notes=[]
 for cid in (selected.combo_id.unique() if len(selected) else []):
  dest=OUT/'review'/cid;dest.mkdir(exist_ok=True);folder=OUT/'combinations'/cid;K=json.loads((folder/'selection.json').read_text())['K'];chosen=selected[selected.combo_id.eq(cid)];units=pd.read_parquet(folder/'spotted_units.parquet');sections=pd.concat([pd.read_parquet(f) for f in (folder/'spotted_sections').glob('*.parquet')],ignore_index=True)
  # Fixed-label leave-one-dataset-out: filter sections THEN reaggregate shared units.
  loo=[];matched=[];dist=[];influence=[]
  for hit in chosen.itertuples():
   ss=sections[(sections.layer==hit.layer)&(sections.track=='SOURCE_P05')&(sections.target_kind=='NICHE')&(sections.target_id==hit.target_id)&(sections.metric==hit.metric)];good=ss[ss.eligible&ss.effect.notna()];ud=good.groupby('analysis_unit_uid').effect.mean();loo.append(dict(target_id=hit.target_id,layer=hit.layer,metric=hit.metric,removed_dataset='NONE',n_units=len(ud),mean_effect=ud.mean(),median_effect=ud.median()))
   for ds in sorted(sections.dataset_id.unique()):
    v=good[good.dataset_id.ne(ds)].groupby('analysis_unit_uid').effect.mean();loo.append(dict(target_id=hit.target_id,layer=hit.layer,metric=hit.metric,removed_dataset=ds,n_units=len(v),mean_effect=v.mean(),median_effect=v.median(),same_direction=bool(v.mean()*hit.mean_effect>0)))
   for uid in ud.index:
    v=ud.drop(uid);influence.append(dict(target_id=hit.target_id,layer=hit.layer,metric=hit.metric,removed_unit=uid,n_remaining=len(v),mean_remaining=v.mean(),same_direction=bool(v.mean()*hit.mean_effect>0),absolute_mean_change=abs(v.mean()-ud.mean())))
   pair=sections[(sections.track=='SOURCE_P05')&(sections.target_kind=='NICHE')&(sections.target_id==hit.target_id)&(sections.metric==hit.metric)&sections.eligible];a=pair[pair.layer==LAYERS[0]];b=pair[pair.layer==LAYERS[1]];joined=a.merge(b,on=['section_id','analysis_unit_uid'],suffixes=('_spot','_epi'),validate='one_to_one');mu=joined.groupby('analysis_unit_uid')[['effect_spot','effect_epi']].mean().reset_index();mu['target_id']=hit.target_id;mu['metric']=hit.metric;matched.append(mu)
   ds=good.groupby('dataset_id').effect.agg(['count','mean','median']);fig,ax=plt.subplots(figsize=(7,4));scopes=['ALL_DATASETS']+sorted(sections.dataset_id.unique());values=[]
   for scope in scopes:
    v=good if scope=='ALL_DATASETS' else good[good.dataset_id.eq(scope)];v=v.groupby('analysis_unit_uid').effect.mean();values.append(v.to_numpy());dist.extend([dict(target_id=hit.target_id,layer=hit.layer,metric=hit.metric,scope=scope,analysis_unit_uid=uid,effect=val) for uid,val in v.items()])
   for i,v in enumerate(values):
    if len(v):ax.scatter(v,np.full(len(v),i),s=12,alpha=.65);ax.plot([np.median(v)]*2,[i-.2,i+.2],color='black')
   ax.axvline(0,color='gray',lw=1);ax.set_yticks(range(len(scopes)),scopes,fontsize=7);ax.set(xlabel='Overlap proportion difference' if hit.metric=='O' else 'Distance difference (spot pitches)',title=f'{cid} {hit.target_id} {hit.metric} {hit.layer}\n{level}; q_global={hit.q_global_BH:.3g}; descriptive selected effects');save(fig,dest/f'{hit.target_id}_{hit.metric}_{hit.layer}_effects')
  table(dest/'leave_one_dataset_out.tsv',loo);table(dest/'leave_one_unit_out.tsv',influence);table(dest/'effect_distribution.parquet',dist);table(dest/'matched_score_units.parquet',pd.concat(matched,ignore_index=True).drop_duplicates() if matched else pd.DataFrame())
  # Show overlap, distance and non-overlap distance separately for each selected target/layer.
  for target,layer in chosen[['target_id','layer']].drop_duplicates().itertuples(index=False,name=None):
   fig,axs=plt.subplots(1,3,figsize=(15,4));scopes=['ALL_DATASETS']+sorted(sections.dataset_id.unique())
   for ax,metric in zip(axs,['O','D','D_out']):
    yy=[]
    for i,scope in enumerate(scopes):
     use=sections[(sections.target_id==target)&(sections.target_kind=='NICHE')&(sections.layer==layer)&(sections.track=='SOURCE_P05')&(sections.metric==metric)&sections.eligible]
     if scope!='ALL_DATASETS':use=use[use.dataset_id.eq(scope)]
     v=use.groupby('analysis_unit_uid').effect.mean();yy.append(scope+' N='+str(len(v)))
     if len(v):ax.scatter(v,np.full(len(v),i),s=9,alpha=.6);ax.plot([v.median()]*2,[i-.2,i+.2],color='black')
    ax.axvline(0,color='gray',lw=.7);ax.set_yticks(range(len(scopes)),yy,fontsize=6);ax.set_title(metric);ax.set_xlabel('Proportion difference' if metric=='O' else 'Spot-pitch difference')
   fig.suptitle(cid+' '+target+' '+layer+'; selected descriptive distributions');fig.tight_layout();save(fig,dest/f'{target}_{layer}_O_D_Dout')
  # Composition heatmap and contributions for all niches.
  prof=pd.read_parquet(folder/'profiles.parquet');mat=prof.groupby(['niche','celltype'])['mean'].mean().unstack().reindex(columns=TYPES);fig,ax=plt.subplots(figsize=(13,4));im=ax.imshow(mat,aspect='auto',cmap='viridis');ax.set_xticks(range(30),TYPES,rotation=90,fontsize=7);ax.set_yticks(range(K),['niche'+str(i) for i in range(K)]);fig.colorbar(im,ax=ax,label='Mean R1/R2 relative composition');save(fig,dest/'composition')
  cov=pd.read_parquet(folder/'coverage.parquet');con=cov.groupby(['niche','dataset_id']).n_accepted.sum().unstack(fill_value=0);fig,ax=plt.subplots(figsize=(8,4));con.div(con.sum(1),axis=0).plot.bar(stacked=True,ax=ax);ax.legend(fontsize=6,bbox_to_anchor=(1.01,1));ax.set(ylabel='Accepted-anchor share');save(fig,dest/'dataset_contribution');table(dest/'unit_contribution.tsv',cov.groupby(['niche','analysis_unit_uid']).n_accepted.sum().reset_index())
  uc=cov.groupby(['niche','analysis_unit_uid']).n_accepted.sum().unstack(fill_value=0);fig,ax=plt.subplots(figsize=(13,4));im=ax.imshow(uc.div(uc.sum(1),axis=0),aspect='auto',cmap='viridis');ax.set_xticks(range(len(uc.columns)),uc.columns,rotation=90,fontsize=4);ax.set_yticks(range(K),['niche'+str(i) for i in range(K)]);ax.set_title('Analysis-unit contributions; not all identities verified');fig.colorbar(im,ax=ax,label='Accepted-anchor share');save(fig,dest/'unit_contribution')
  if not (dest/'reference_embedding.npz').exists():subprocess.run(['/usr/bin/python3',str(OUT/'code/refit.py'),cid,'--rebuild'],check=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
  if not (dest/'umap.npz').exists():subprocess.run(['/usr/bin/python3',str(OUT/'code/umap_display.py'),cid],check=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
  um=np.load(dest/'umap.npz');fig,ax=plt.subplots(figsize=(6,5));ax.scatter(um['xy'][:,0],um['xy'][:,1],c=um['labels'],s=1,cmap='tab20');ax.set(title=cid+' reference UMAP (display only)',xlabel='UMAP1',ylabel='UMAP2');save(fig,dest/'reference_umap')
  labels=np.full(len(master),-1,np.int16);ll=np.load(folder/'labels.npz');labels[ll['anchor_index']]=ll['accepted'];scov=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');members=c.set_index('combo_id').loc[cid,'datasets'].split('|');chosen_sections=[];sens=[];display=[]
  for ds in members:
   secds=scov[scov.dataset_id.eq(ds)].section_id.tolist()
   if ds==INHOUSE:chosen_sections+=secds
   else:
    covs=sections[sections.dataset_id.eq(ds)].groupby('section_id').n_source_supported.max();med=covs.median();chosen_sections+=sorted(secds,key=lambda s:(-int(covs.get(s,0)>=30),abs(covs.get(s,0)-med),s))[:2]
  for sec in scov[scov.dataset_id.isin(members)].section_id:
   t,a,g=section(sec,labels,master,base)
   for k in [6,10,18]:
    src,info=source(sec,t,a,g,k,folder);se=effects(sec,t,a,g,src,K,k,domains=['MAIN','GEOMETRY_LOCAL_SUPPORT'] if k==10 else ['MAIN']);keys={(r.target_id,r.layer,r.metric) for r in chosen.itertuples()};se=se[[ (r.target_id,r.layer,r.metric) in keys for r in se.itertuples() ]];sens.append(se)
   if sec in chosen_sections:
    for hit in chosen.itertuples():
     fig,axs=plt.subplots(1,3,figsize=(10,3));xy=a['xy_pitch'];target=t.accepted_label.eq(int(hit.target_id.replace('niche',''))).to_numpy();src,info=source(sec,t,a,g,10,folder);li=LAYERS.index(hit.layer);hot=src['node_index'][src[f'H{li}_0']]
     for ax in axs:ax.scatter(xy[:,0],xy[:,1],s=1,c='#dddddd');ax.set_aspect('equal');ax.invert_yaxis();ax.set_axis_off()
     axs[0].scatter(xy[target,0],xy[target,1],s=2,c='#377eb8');axs[0].set_title(hit.target_id)
     valid=np.isfinite(t[hit.layer]);axs[1].scatter(xy[valid,0],xy[valid,1],c=t.loc[valid,hit.layer],s=2,cmap='viridis');axs[1].set_title('MP10 existing score')
     axs[2].scatter(xy[hot,0],xy[hot,1],s=2,c='#e41a1c');axs[2].set_title(f'Formal source hot: {len(hot)}')
     fig.suptitle(sec);save(fig,dest/f'spatial_{sec}_{hit.target_id}_{hit.metric}_{hit.layer}')
    display.append(dict(section_id=sec,selection='all inhouse' if sec.startswith('Inhouse') else 'coverage >=30 then nearest median source support, section ID tie',n_E0=int(t.E0.sum()),n_accepted=int(t.accepted_label.ge(0).sum())))
  table(dest/'section_display_selection.tsv',display);table(dest/'spotted_sensitivity.parquet',pd.concat(sens,ignore_index=True));ringrun(cid,sets=['COMPLETE_RINGS','TME_MATCHED_RAW','TME_CONDITIONAL','SPARSE_BLOCK'])
  # Explicit matched-composition comparison against frozen-label controls.
  mapping=pd.read_csv(folder/'old_niche_mapping.tsv',sep='\t');frozen=pd.read_parquet(folder/'frozen_tests.parquet');compare=[]
  for hit in chosen.itertuples():
   k=int(hit.target_id.replace('niche',''));mp=mapping[mapping.niche.eq(k)].sort_values(['profile_L2','n_overlap'],ascending=[True,False]).iloc[0];old=frozen[(frozen.scope=='ALL_DATASETS')&(frozen.track=='SOURCE_P05')&(frozen.target_kind=='NICHE')&(frozen.target_id==f'niche{int(mp.old_niche)}')&(frozen.metric==hit.metric)&(frozen.layer==hit.layer)&(frozen.subset=='ALL_POLICY_UNITS')].iloc[0];compare.append(dict(target=hit.target_id,old_target=old.target_id,profile_L2=mp.profile_L2,overlap_fraction=mp.overlap_fraction,layer=hit.layer,metric=hit.metric,refit_effect=hit.mean_effect,frozen_effect=old.mean_effect,refit_N=hit.n_units,frozen_N=old.n_units,same_direction=hit.mean_effect*old.mean_effect>0,warning='Composition nearest mapping may be split/merge; old source E0 differs; not independent validation'))
  table(dest/'refit_frozen_comparison.tsv',compare);notes.append(dict(combo_id=cid,n_relations=len(chosen),figures=len(list(dest.glob('*.png'))),bytes=bytes_used(dest)))
 table(OUT/'review/review_manifest.tsv',notes);js(OUT/'review/COMPLETE.json',dict(status='COMPLETE',time=time.time(),level=level,selected_combinations=len(notes)))
if __name__=='__main__':run()
