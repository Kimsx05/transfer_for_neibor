from pathlib import Path
import numpy as np,pandas as pd,json,zipfile,importlib.metadata as md
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score,adjusted_mutual_info_score
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[1];sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];types=json.loads((R/'type_mapping.json').read_text());F=R/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':12,'lines.linewidth':1.1,'pdf.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white'})
metrics=[];comp=[];quality=[];coverage=[];labels={};native=[];validation=[];summaries=[]
def export(fig,stem,handles=None):
 fig.savefig(F/(stem+'.png'),dpi=300,bbox_inches='tight');fig.savefig(F/(stem+'.pdf'),bbox_inches='tight')
 if handles:
  for ax in fig.axes:
   leg=ax.get_legend()
   if leg:leg.remove()
  for leg in list(fig.legends):leg.remove()
  fig.savefig(F/(stem+'_panel_only.pdf'),bbox_inches='tight');
  from matplotlib.backends.backend_pdf import PdfPages
  with PdfPages(F/(stem+'_legend_only.pdf')) as pdf:
   for start in range(0,len(handles),100):
    lf,la=plt.subplots(figsize=(12,10));la.axis('off');la.legend(handles=handles[start:start+100],loc='center',frameon=False,ncol=4,fontsize=9);pdf.savefig(lf,bbox_inches='tight');plt.close(lf)
 plt.close(fig)
for ref,cats in types.items():
 for radius in [50,200]:
  o=R/f'{ref}_r{radius}';j=json.loads((o/'complete.json').read_text());summaries.append(j);d=pd.read_parquet(o/'labels.parquet');labels[ref,radius]=d
  p=pd.concat([pd.read_parquet(R/s/f'r{radius}.parquet') for s in sections],ignore_index=True);assert d[['section','bin_id']].equals(p[['section','bin_id']]);assert not d.duplicated(['section','bin_id']).any()
  tm=pd.read_csv(o/'training_anchor_ids.tsv',sep='\t');keys=pd.MultiIndex.from_frame(d[['section','bin_id']]);train=keys.isin(pd.MultiIndex.from_frame(tm));assert (d.loc[train,'A_assignment_source']=='native_train').all() and (d.loc[train,'B_assignment_source']=='native_train').all()
  for scope,mask in [('native_T',train),('full_S',d.evaluable.to_numpy()),('full_assigned_both',d.evaluable.to_numpy()&(d.A_niche!='unresolved')&(d.B_niche!='unresolved'))]:
   a=d.loc[mask,'A_niche'];b=d.loc[mask,'B_niche'];metrics.append(dict(reference=ref,radius=radius,scope=scope,n=len(a),ARI=adjusted_rand_score(a,b),NMI=normalized_mutual_info_score(a,b),AMI=adjusted_mutual_info_score(a,b)))
   pd.crosstab(a,b).to_csv(o/f'method_crosstab_{scope}.tsv',sep='\t')
  coverage.append(pd.read_csv(o/'coverage.tsv',sep='\t').assign(reference=ref,radius=radius))
  raw=pd.concat([pd.read_parquet(R/s/f'neighbor_marker_r{radius}.parquet') for s in sections],ignore_index=True);assert raw.bin_id.tolist()==d.bin_id.tolist()
  for method in ['A','B']:
   joined=pd.concat([d,p.drop(columns=['section','bin_id']),raw.drop(columns='bin_id')],axis=1)
   for cluster,g in joined.groupby(method+'_niche'):
    cs=g[[ref+'_P_'+t for t in cats]].mean();row=dict(reference=ref,radius=radius,method=method,niche=cluster,n=len(g));row.update({t:cs[ref+'_P_'+t] for t in cats});comp.append(row)
    quality.append(dict(reference=ref,radius=radius,method=method,niche=cluster,n=len(g),known_median=g[ref+'_n_known'].median(),known_fraction_median=g[ref+'_known_fraction'].median(),coverage_q10=g.input_coverage_component.quantile(.1),low_information_mean=g.low_information_fraction.mean(),reference_conflict_mean=g.reference_conflict_fraction.mean(),**{c:g[c].mean() for c in raw.columns if c!='bin_id'}))
   c=pd.DataFrame([x for x in comp if x['reference']==ref and x['radius']==radius and x['method']==method]).set_index('niche');c=c.sort_values('n',ascending=False)
   if len(c)>32:
    tail=c.iloc[31:];other={t:np.average(tail[t],weights=tail['n']) for t in cats};other['n']=tail['n'].sum();c=pd.concat([c.iloc[:31],pd.DataFrame([other],index=['Other IDs (display only)'])])
   fig,ax=plt.subplots(figsize=(max(6,len(cats)*.8),max(3,len(c)*.24)));im=ax.imshow(c[cats],vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set_xticks(range(len(cats)),cats,rotation=45,ha='right');ax.set_yticks(range(len(c)),[method+' '+x for x in c.index]);fig.colorbar(im,ax=ax,label='Mean known-neighbor fraction');export(fig,f'composition_{ref}_{radius}_{method}')
  validation.append(dict(reference=ref,radius=radius,unique_IDs=True,common_P_S_T=True,API_verified=j['API_composition_verified'],native_not_overwritten=True,all_shared_anchors=len(d)))
for name,data in [('method_comparison',metrics),('niche_composition',comp),('niche_quality_marker_support',quality),('selected_models',summaries),('validation',validation)]:pd.DataFrame(data).to_csv(R/'tables'/f'{name}.tsv',sep='\t',index=False)
pd.concat(coverage).to_csv(R/'tables/coverage.tsv',sep='\t',index=False)
# Paired scale/reference comparisons retain common evaluable denominator explicitly.
paired=[]
for ref in types:
 a,b=labels[ref,50],labels[ref,200];assert a[['section','bin_id']].equals(b[['section','bin_id']]);ok=a.evaluable&b.evaluable
 for method in ['A','B']:
  x,y=a.loc[ok,method+'_niche'],b.loc[ok,method+'_niche'];pd.crosstab(x,y).to_csv(R/'tables'/f'scale_{ref}_{method}.tsv',sep='\t');paired.append(dict(comparison='scale',reference=ref,method=method,n=int(ok.sum()),ARI=adjusted_rand_score(x,y),NMI=normalized_mutual_info_score(x,y)))
for radius in [50,200]:
 a,b=labels['blca',radius],labels['hpca',radius];ok=a.evaluable&b.evaluable
 for method in ['A','B']:
  x,y=a.loc[ok,method+'_niche'],b.loc[ok,method+'_niche'];pd.crosstab(x,y).to_csv(R/'tables'/f'reference_r{radius}_{method}.tsv',sep='\t');paired.append(dict(comparison='reference',reference=str(radius),method=method,n=int(ok.sum()),ARI=adjusted_rand_score(x,y),NMI=normalized_mutual_info_score(x,y)))
pd.DataFrame(paired).to_csv(R/'tables/paired_comparisons.tsv',sep='\t',index=False)
for s in sections:
 nodes=pd.read_parquet(R/s/'nodes.parquet');a=nodes[nodes.epi_anchor]
 for ref in types:
  fig,axs=plt.subplots(2,2,figsize=(12,11));handles=[]
  for row,method in enumerate(['A','B']):
   for col,radius in enumerate([50,200]):
    ax=axs[row,col];z=labels[ref,radius];z=z[z.section==s].set_index('bin_id').loc[a.bin_id];labs=sorted(set(labels[ref,radius][method+'_niche'])-{'unresolved'});pal=plt.get_cmap('hsv');colors={v:pal(((i*.61803398875)+(0 if method=='A' else .23+radius/1000+(0.11 if ref=='hpca' else 0)))%1) for i,v in enumerate(labs)};
    if method=='A' and len(labs)==2:
     cc=[v for v in comp if v['reference']==ref and v['radius']==radius and v['method']=='A' and v['niche']!='unresolved'];high=max(cc,key=lambda v:v['Epithelial'])['niche'];colors={k:('#f5b264' if k==high else '#8752a3') for k in labs}
    colors['unresolved']='#444444'
    ax.scatter(nodes.x_um,nodes.y_um,s=.12,c='#e8e0d0',rasterized=True);ax.scatter(a.x_um,a.y_um,c=[colors[x] for x in z[method+'_niche']],s=.4,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel(f'{method}, r={radius} μm: x (μm)');ax.set_ylabel('Relative y (μm)' if s.startswith('MPBC') else 'y (μm)')
    hh=[Line2D([],[],marker='o',linestyle='',color=colors[k],label=f'{method} r{radius}: {k}') for k in labs];handles+=hh
  handles+=[Line2D([],[],marker='o',linestyle='',color='#444444',label='Unresolved anchor'),Line2D([],[],marker='o',linestyle='',color='#e8e0d0',label='Supplied background (not niche)')];fig.legend(handles=[h for prefix in ['A r50:', 'A r200:', 'B r50:', 'B r200:'] for h in [v for v in handles if v.get_label().startswith(prefix)][:6]]+handles[-2:],loc='center left',bbox_to_anchor=(1,.5),fontsize=8,title='Up to6 keys/panel; complete legend in PDF');export(fig,f'{s}_{ref}_methods_scales',handles)
versions={}
for pkg in ['squidpy','scanpy','anndata','scipy','scikit-learn','igraph','leidenalg','pandas','numpy']:
 try:versions[pkg]=md.version(pkg)
 except:versions[pkg]='unavailable'
(R/'versions.json').write_text(json.dumps(versions,indent=2));print('FINAL_TABLES_FIGURES_COMPLETE',flush=True)
