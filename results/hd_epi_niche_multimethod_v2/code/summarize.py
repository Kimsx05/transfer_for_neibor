from common import *
from sklearn.metrics import adjusted_rand_score as ari,normalized_mutual_info_score as nmi
from itertools import combinations

def eta(y,v):
 good=(y>=0)&np.isfinite(v);y=y[good];v=v[good]
 if len(v)<2 or np.var(v)==0:return np.nan
 tot=np.sum((v-v.mean())**2);return float(sum(len(v[y==k])*(v[y==k].mean()-v.mean())**2 for k in np.unique(y))/tot)
if __name__=='__main__':
 a=anchors();T=train_index(a);labels={};selected=[];cand=[];coverage=[];comps=[];quality=[];technical=[];diag=[];ringsum=[];comparison=[];index=[]
 root=R/'all_anchor_assignments';root.mkdir(exist_ok=True)
 ds={r:pd.concat([pd.read_parquet(V/s/f'r{r}.parquet') for s in SECTIONS],ignore_index=True) for r in [50,200]}
 raw={r:pd.concat([pd.read_parquet(V/s/f'neighbor_marker_r{r}.parquet') for s in SECTIONS],ignore_index=True).drop(columns='bin_id') for r in [50,200]}
 for r in [50,200]:
  extra=pd.concat([pd.read_parquet(R/'tables'/f'{s}_raw_panel_neighbors_r{r}.parquet') for s in SECTIONS],ignore_index=True).drop(columns=['bin_id','section']);raw[r]=pd.concat([raw[r],extra],axis=1)
 for o in sorted((R/'models').iterdir()):
  if not o.is_dir() or not (o/'assignments.parquet').exists():continue
  z=pd.read_parquet(o/'assignments.parquet');assert np.array_equal(a.key,z.section+'::'+z.bin_id);j=json.loads((o/'complete.json').read_text());j['tag']=o.name;selected.append(j);labels[o.name]=z;index.append(dict(tag=o.name,path=str((o/'assignments.parquet').relative_to(R))))
  link=root/(o.name+'.parquet')
  if not link.exists():link.symlink_to(Path('../models')/o.name/'assignments.parquet')
  c=pd.read_csv(o/'candidate_metrics.tsv',sep='\t');c['tag']=o.name;cand.append(c)
  for (s,src),g in z.groupby(['section','source']):coverage.append(dict(tag=o.name,section=s,source=src,n=len(g),anchors=int((a.section==s).sum()),evaluable=int(z.loc[z.section==s,'evaluable'].sum())))
  r=j['radius'];refs=TYPES if j['reference']=='both' else [j['reference']];y=z.label.to_numpy();d=ds[r]
  for ref in refs:
   cp=d[[ref+'_P_'+t for t in TYPES[ref]]].copy();cp.columns=TYPES[ref];cp['label']=y
   for lab,g in cp.groupby('label'):
    rec=dict(tag=o.name,reference=ref,radius=r,label=int(lab),n=len(g));rec.update(g.drop(columns='label').mean().to_dict());rec.update({k+'_sd':float(v) for k,v in g.drop(columns='label').std().items()});comps.append(rec)
  q=d[['nCount_RNA','nFeature_RNA','n_analysis_neighbors','blca_n_known','hpca_n_known','low_information_fraction','reference_conflict_fraction','input_coverage_component','blca_known_fraction','hpca_known_fraction']].copy();q=pd.concat([q,raw[r]],axis=1);q['label']=y
  for lab,g in q.groupby('label'):
   row=dict(tag=o.name,label=int(lab),n=len(g));row.update(g.drop(columns='label').mean().to_dict());quality.append(row)
  tt=dict(tag=o.name,section_NMI=float(nmi(a.section[y>=0],y[y>=0])))
  for col in q.columns.drop('label'):tt['eta2_'+col]=eta(y,q[col].to_numpy())
  technical.append(tt)
  gd=o/'graph_diagnostics.tsv'
  if gd.exists():gg=pd.read_csv(gd,sep='\t');gg['tag']=o.name;diag.append(gg)
  if j['method']=='MENDER_adapted':
   ref=j['reference'];cs=[];qc=[]
   for s in SECTIONS:
    f=np.load(R/'models'/f'{s}_rings.npz');ix=[list(f['columns']).index(ref+'::'+t) for t in TYPES[ref]];cs.append(f['counts'][:,:r//25,ix]);qc.append(f['counts'][:,:r//25,1])
   cs=np.concatenate(cs);qc=np.concatenate(qc)
   for lab in np.unique(y):
    mask=y==lab
    for ring in range(r//25):
     known=cs[mask,ring].sum(1);den=qc[mask,ring];props=np.divide(cs[mask,ring],known[:,None],out=np.full(cs[mask,ring].shape,np.nan),where=known[:,None]>0);row=dict(tag=o.name,label=int(lab),ring=ring,inner_um=ring*25,outer_um=(ring+1)*25,mean_QC=float(den.mean()),mean_known=float(known.mean()),mean_Unknown=float((den-known).mean()),empty_known_fraction=float((known==0).mean()));row.update(dict(zip(TYPES[ref],np.nanmean(props,axis=0))));ringsum.append(row)
 def pair(aa,bb,kind):
  if aa not in labels or bb not in labels:return
  x,y=labels[aa],labels[bb]
  for scope in ['common_native_T','common_evaluable','common_assigned']:
   ix=np.zeros(len(a),bool);ix[T]=True
   if scope!='common_native_T':ix=x.evaluable.to_numpy()&y.evaluable.to_numpy()
   if scope=='common_assigned':ix&=(x.label>=0)&(y.label>=0)
   comparison.append(dict(kind=kind,left=aa,right=bb,scope=scope,n=int(ix.sum()),fraction_shared_anchors=float(ix.mean()),ARI=ari(x.loc[ix,'label'],y.loc[ix,'label']),NMI=nmi(x.loc[ix,'label'],y.loc[ix,'label'])))
 for ref in TYPES:
  for r in [50,200]:
   tags=[f'hellinger_kmeans_{ref}_r{r}_NONE',f'adapted_leiden_{ref}_r{r}_NONE',f'MENDER_adapted_{ref}_r{r}_ring_counts',f'CellCharter_both_r{r}_scVI']
   for aa,bb in combinations(tags,2):pair(aa,bb,'method')
   for m in ['hellinger_kmeans','adapted_leiden']:pair(f'{m}_{ref}_r{r}_NONE',f'{m}_{ref}_r{r}_WEAK','Harmony')
  for m,rep in [('hellinger_kmeans','NONE'),('adapted_leiden','NONE'),('MENDER_adapted','ring_counts')]:pair(f'{m}_{ref}_r50_{rep}',f'{m}_{ref}_r200_{rep}','scale')
 for r in [50,200]:
  for m,rep in [('hellinger_kmeans','NONE'),('adapted_leiden','NONE'),('MENDER_adapted','ring_counts')]:pair(f'{m}_blca_r{r}_{rep}',f'{m}_hpca_r{r}_{rep}','reference')
  pair(f'CellCharter_both_r{r}_scVI',f'CellCharter_both_r{r}_scVI_without_Top20','feature_sensitivity')
 pair('CellCharter_both_r50_scVI','CellCharter_both_r200_scVI','scale')
 pair('MENDER_adapted_blca_r50_ring_counts','MENDER_adapted_blca_r50_per_ring_normalized','ring_support_sensitivity')
 for name,data in [('selected_models',selected),('coverage',coverage),('niche_composition',comps),('quality_marker_support',quality),('technical_association',technical),('MENDER_ring_composition',ringsum),('method_comparisons',comparison),('assignment_index',index)]:pd.DataFrame(data).to_csv(R/'tables'/f'{name}.tsv',sep='\t',index=False)
 if cand:pd.concat(cand).to_csv(R/'tables/candidate_metrics.tsv',sep='\t',index=False)
 if diag:pd.concat(diag).to_csv(R/'tables/identical_input_graph_diagnostics.tsv',sep='\t',index=False)
 scopes=[]
 for j in selected:scopes.append(dict(tag=j['tag'],niche_fit='T_v2 first four sections only',MPBC2_in_niche_fit=False,MPBC2_representation=('all_five_sections_transductive Harmony' if j['representation']=='WEAK' else 'official scArches query adaptation' if j['method']=='CellCharter' else 'not used in representation fitting; transform only')))
 pd.DataFrame(scopes).to_csv(R/'tables/integration_scope.tsv',sep='\t',index=False);print('SUMMARY_COMPLETE',len(selected),flush=True)

# Explain finer labels at exactly shared cumulative P, without choosing centers/models from outcomes.
extra=[];cases=[]
for ref in TYPES:
 for r in [50,200]:
  d,P,ok=composition(ref,r);u,inv,frequency=np.unique(P[T],axis=0,return_inverse=True,return_counts=True);marker=raw[r].iloc[T].reset_index(drop=True)
  for tag,z in labels.items():
   if int(z.radius.iloc[0])!=r or str(z.reference.iloc[0]) not in [ref,'both']:continue
   y=z.iloc[T].label.to_numpy();pairs=pd.DataFrame({'profile':inv,'label':y});ng=pairs.groupby('profile').label.nunique();multi=ng[ng>1].index.to_numpy();extra.append(dict(tag=tag,composition_reference=ref,exact_P_profiles=len(u),profiles_split_across_working_labels=len(multi),training_rows_in_split_profiles=int(frequency[multi].sum()),note='Additional representation/correction/partition distinction; requires raw/radial evidence, not automatically biological'))
   eligible=[i for i in multi if frequency[i]>=50];eligible=sorted(eligible,key=lambda i:-frequency[i])[:5]
   for profile in eligible:
    ids=np.flatnonzero(inv==profile)
    for lab in np.unique(y[ids]):
     ix=ids[y[ids]==lab];row=dict(tag=tag,composition_reference=ref,profile=int(profile),profile_frequency=int(frequency[profile]),label=int(lab),n=len(ix));row.update(dict(zip(TYPES[ref],u[profile])));row.update(marker.iloc[ix].mean().to_dict());row['section_counts']=str(a.iloc[T[ix]].section.value_counts().to_dict());cases.append(row)
pd.DataFrame(extra).to_csv(R/'tables/same_cumulative_P_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(cases).to_csv(R/'tables/same_P_marker_case_summaries.tsv',sep='\t',index=False)

from extended_diagnostics import run as extended_diagnostics
extended_diagnostics()

import runpy
runpy.run_path(str(R/'code/diagnose_ring_support.py'),run_name='__main__')
