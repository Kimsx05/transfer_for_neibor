exec(open('mp10_subtype_screen_round1/scripts/05_abundance_figures.py').read().split('comps=')[0])
e=pd.read_csv(T/'CellChat_detected_edges.tsv.gz',sep='\t');a=pd.read_csv(T/'CellChat_major_aggregated.tsv',sep='\t');cv=pd.read_csv(T/'CellChat_group_coverage.tsv',sep='\t')
map_=dict(pd.read_csv(T/'subtype_parent.tsv',sep='\t')[['subtype','major']].values);map_.update({x:'Epithelial' for x in ['Malignant_like_epithelial','Normal_epithelial','Other_epithelial_unresolved']})
cv['major']=cv.communication_type.map(map_);e['source_major']=e.source.map(map_);e['target_major']=e.target.map(map_)
majors=sorted(set(map_.values()));io=[];pws=[]
specs={'Normal_Tumor':['Normal','Tumor'],'T1_T2':['T1','T2'],'MP10_low_high':['low','high']}
for (ds,co),z in a.groupby(['scope','comparison']):
 gs=specs[co];mats={}
 for metric,title in [('total_strength','Sum of subtype-LR strengths'),('interaction_count','Number of detected subtype-LR edges')]:
  for g in gs:
   mat=z[z.condition.eq(g)].pivot(index='source_major',columns='target_major',values=metric).reindex(index=majors,columns=majors).fillna(0)
   valid=set(cv[cv.scope.eq(ds)&cv.comparison.eq(co)&cv.condition.eq(g)&cv.evaluable].major)
   for label in set(majors)-valid:mat.loc[label,:]=np.nan;mat.loc[:,label]=np.nan
   mats[g]=mat
  vals=np.concatenate([q.values.ravel() for q in mats.values()]);vmax=np.nanmax(vals) if np.isfinite(vals).any() else 1
  diff=mats[gs[1]]-mats[gs[0]];dmax=np.nanmax(np.abs(diff.values)) if np.isfinite(diff.values).any() else 1
  fig,axs=plt.subplots(1,3,figsize=(15,5.4))
  for j,(g,mat) in enumerate(list(mats.items())+[('Difference',diff)]):
   im=axs[j].imshow(np.ma.masked_invalid(mat.values),cmap='RdBu_r' if j==2 else 'viridis',vmin=-max(dmax,1e-10) if j==2 else 0,vmax=max(dmax if j==2 else vmax,1e-10));axs[j].set_xticks(range(len(majors)),majors,rotation=60,ha='right',fontsize=8);axs[j].set_yticks(range(len(majors)),majors,fontsize=8);axs[j].set(title=g,xlabel='Receiver',ylabel='Sender' if j==0 else '');fig.colorbar(im,ax=axs[j],shrink=.65)
   if not np.isfinite(mat.values).any():axs[j].text(.5,.5,'Not evaluable',transform=axs[j].transAxes,ha='center',va='center',fontsize=11)
   axs[j].text(.5,-.35,'Blank cells: not evaluable',transform=axs[j].transAxes,ha='center',fontsize=8)
  fig.suptitle(f'{ds}: {co}\n{title} — aggregated from subtype networks',y=1.02);fig.tight_layout();save(fig,ds,'network_major_'+co+'_'+metric)
 for g,q in z.groupby('condition'):
  for typ in majors:io.append(dict(scope=ds,comparison=co,condition=g,major=typ,outgoing=q.loc[q.source_major.eq(typ),'total_strength'].sum(),incoming=q.loc[q.target_major.eq(typ),'total_strength'].sum()))
# Explicit signalling category and direction retained for every pathway, not just displayed top pathways.
y=e[e.source_major.eq('Epithelial')|e.target_major.eq('Epithelial')].copy();y['direction']=np.where(y.source_major.eq('Epithelial'),'Epithelial -> TME','TME -> Epithelial');y=y[y.source_major.ne(y.target_major)]
pw=y.groupby(['scope','comparison','condition','category','direction','pathway']).agg(strength=('probability','sum'),n_detected_edges=('LR_ID','size')).reset_index();pw.to_csv(T/'CellChat_epithelial_pathways.tsv',sep='\t',index=False)
pd.DataFrame(io).to_csv(T/'CellChat_major_incoming_outgoing.tsv',sep='\t',index=False)
for (ds,co),z in pw.groupby(['scope','comparison']):
 order=z.groupby('pathway').strength.sum().nlargest(20).index.tolist();z=z[z.pathway.isin(order)]
 fig,axs=plt.subplots(1,2,figsize=(10,max(5,len(order)*.27)),sharey=True)
 for ax,di in zip(axs,['Epithelial -> TME','TME -> Epithelial']):
  q=z[z.direction.eq(di)];xx=q.condition.map({v:i for i,v in enumerate(specs[co])});yy=q.pathway.map({v:i for i,v in enumerate(order)})
  ax.scatter(xx,yy,s=20+q.n_detected_edges*3,c=q.strength,cmap='viridis',vmin=0,vmax=z.strength.max());ax.set_xticks([0,1],specs[co]);ax.set_yticks(range(len(order)),order);ax.set(title=di,xlim=(-.5,1.5));ax.grid(alpha=.15)
 fig.suptitle(f'{ds}: {co} epithelial pathways\nColor: summed strength; area: detected subtype-LR edges');fig.tight_layout();save(fig,ds,'pathway_bubbles_'+co)
print('Network overviews complete')
