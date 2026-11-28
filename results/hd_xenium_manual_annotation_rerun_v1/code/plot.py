from data import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import colorsys,itertools,math
p=sys.argv[1];d=O/p;out=O/'figures'/p;out.mkdir(exist_ok=True)
z=pd.read_parquet(d/'objects.parquet');qc=z
markers={'T_NK':'PTPRC CD3D CD3E CD3G TRAC CD2 CD4 CD8A CD8B NKG7 KLRD1 GNLY','Treg':'FOXP3 IL2RA CTLA4 IKZF2 TIGIT','B_Plasma':'MS4A1 CD79A CD79B JCHAIN MZB1 SDC1','Myeloid':'LYZ LST1 TYROBP CSF1R CD68 FCGR3A S100A8 S100A9 FOLR2 SPP1 LAMP3','Stromal':'COL1A1 COL1A2 DCN LUM PDGFRA DPT ACTA2 RGS5 PDGFRB MYH11','Endothelial':'PECAM1 VWF CLDN5 KDR EMCN ACKR1 PROX1 LYVE1','Epithelial':'EPCAM KRT8 KRT18 KRT19 KRT5 KRT14 UPK1A UPK2','Mast':'KIT TPSAB1 TPSB2'}
genes=pd.read_csv(d/'features.tsv',sep='\t');ids=genes.feature_id;symbol=genes.gene_symbol
avail=[];by_marker={}
for family,ls in markers.items():
 keep=[]
 for g in ls.split():
  hit=genes[genes.gene_symbol==g];avail.append((p,family,g,'|'.join(hit.feature_id),len(hit),'|'.join(S[p])))
  keep.extend(hit.feature_id)
 by_marker[family]=keep
write(pd.DataFrame(avail,columns=['platform','family','gene_symbol','feature_ids','n_measured_features','samples_assessed']),d/'marker_availability.tsv')
def palette(n):
 first=['#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FDBF6F','#FF7F00','#FB9A99']
 return (first+[plt.matplotlib.colors.to_hex(colorsys.hls_to_rgb((i*.61803398875+.05)%1,.34+.18*(i%3)/2,.65)) for i in range(n) ])[:n]
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=170,bbox_inches='tight')
 plt.close(fig)
def scatter(sub,col,colors,name,title,spatial=False):
 fig,ax=plt.subplots(figsize=(6.5,5.4));key=['x_um','y_um'] if spatial else ['UMAP1','UMAP2']
 a=sub.sample(frac=1,random_state=1729);cc=a[col].map(colors);assert cc.notna().all()
 ax.scatter(a[key[0]],a[key[1]],s=.48,c=cc,linewidths=0,rasterized=True)
 if not spatial:
  for k,v in sub.groupby(col):ax.text(v.UMAP1.median(),v.UMAP2.median(),str(k),ha='center',va='center',fontsize=7,fontweight='bold')
  for axis,k in zip([ax.set_xlim,ax.set_ylim],key):
   lo,hi=z[k].agg(['min','max']);pad=(hi-lo)*.035;axis(lo-pad,hi+pad)
  ax.set_xticks([]);ax.set_yticks([]);ax.spines[:].set_visible(False)
 else:ax.invert_yaxis()
 ax.set_xlabel('x (µm)' if spatial else 'UMAP1');ax.set_ylabel('y (µm)' if spatial else 'UMAP2');ax.set_aspect('equal')
 ax.legend(handles=[Line2D([],[],color=v,ls='',marker='o',ms=4,label=str(k)) for k,v in colors.items()],bbox_to_anchor=(1,.5),loc='center left',ncol=max(1,math.ceil(len(colors)/24)),frameon=False,fontsize=7)
 ax.set_title(title,fontsize=10);save(fig,name)
scatter(z,'sample',dict(zip(S[p],palette(len(S[p])))),'UMAP_sample',p+' · sample')
template=None
for r in R:
 col=cname(r);dr=d/('r'+r);keys=sorted(z[col].unique());colors=dict(zip(keys,palette(len(keys))));write(pd.DataFrame({'cluster':keys,'color':list(colors.values())}),dr/'colors.tsv')
 ct=pd.crosstab(z[col],z['sample']).reindex(columns=S[p],fill_value=0);frac=ct.div(ct.sum(1),axis=0)
 info=pd.DataFrame({'platform':p,'resolution':r,'cluster':ct.index,'n_objects':ct.sum(1).values,'n_samples_present':(ct>0).sum(1).values,'dominant_sample':ct.idxmax(1).values,'dominant_sample_fraction':frac.max(1).values})
 for s in S[p]:info[s+'_n']=ct[s].values;info[s+'_fraction']=frac[s].values
 info['median_GEX_UMI']=z.groupby(col).n_counts.median().values;info['median_nFeature']=z.groupby(col).nFeature.median().values;info['median_area_um2']=z.groupby(col).area_um2.median().values
 info['GEX_secondary_n']=z.source.eq('GEX_secondary').groupby(z[col]).sum().values
 write(info,dr/'cluster_composition.tsv');scatter(z,col,colors,'UMAP_r'+r,p+' · Leiden '+r)
 if r=='.5' or r=='0.5':
  for s in S[p]:
   sub=z[z['sample']==s];scatter(sub,col,colors,s+'_UMAP_r0.5',s+' · joint UMAP');scatter(sub,col,colors,s+'_spatial_r0.5',s+' · Leiden 0.5',True)
  template=info.copy()
  for direction in ['positive','negative']:
   top=pd.read_csv(dr/('top20_'+direction+'.tsv'),sep='\t')
   for n in [10,20]:
    summary=top[top['rank']<=n].groupby('cluster').gene_symbol.agg(', '.join)
    template[f'Top{n}_'+direction]=template.cluster.map(summary)
  for k in ['manual_major','manual_subtype','manual_label','supporting_markers','contradictory_markers','confidence','reviewer','notes']:template[k]=''
  write(template,d/'manual_review.tsv')
cols=['cluster','feature_id','gene_symbol','mean_norm_in','pct_in']
pos=pd.read_csv(d/'r0.5/top20_positive.tsv',sep='\t')
# Prioritize three positive features per cluster for the expression panels, preserving the full Top20 table separately.
top=list(dict.fromkeys(pos.loc[pos['rank']<=3,'feature_id']));want=set(top+sum(by_marker.values(),[]));frames=[]
for b in pd.read_csv(d/'r0.5/pooled_all_genes.tsv.gz',sep='\t',usecols=cols,chunksize=60000):frames.append(b[b.feature_id.isin(want)])
small=pd.concat(frames);write(small,d/'dotplot_values.tsv.gz')
def matrixplot(ls,name,heat=False):
 if not ls:return
 ids_to_symbol=dict(zip(ids,symbol));dup=set(symbol[symbol.duplicated(False)]);cl=sorted(z.leiden_r0_5.unique())
 for j in range(0,len(ls),22):
  page=ls[j:j+22];qq=small[small.feature_id.isin(page)]
  for h in range(0,len(cl),36):
   rows=cl[h:h+36];m=qq.pivot(index='cluster',columns='feature_id',values='mean_norm_in').loc[cl,page];full=np.log1p(m.values);v=full[h:h+len(rows)];det=qq.pivot(index='cluster',columns='feature_id',values='pct_in').loc[rows,page].values
   fig,ax=plt.subplots(figsize=(max(5,.42*len(page)+1.7),max(3,.25*len(rows)+1.8)))
   if heat:
    v=(v-full.mean(axis=0))/np.maximum(full.std(axis=0),1e-8);art=ax.imshow(v,aspect='auto',vmin=-2.5,vmax=2.5,cmap='RdBu_r');legend='gene-wise z of log1p(mean CP10k), all clusters'
   else:
    cy,cx=np.indices(v.shape);art=ax.scatter(cx.ravel(),cy.ravel(),s=det.ravel()*90,c=v.ravel(),cmap='viridis',vmin=0,vmax=max(.1,np.log1p(small.mean_norm_in.max())));ax.invert_yaxis();legend='log1p(mean CP10k)'
    handles=[ax.scatter([],[],s=90*f,color='gray',label=f'{f:.0%}') for f in [.1,.5,1.]];ax.legend(handles=handles,title='Detected',loc='upper left',bbox_to_anchor=(1.22,1),fontsize=7,frameon=False)
   fig.colorbar(art,ax=ax,shrink=.7,label=legend);ax.set_yticks(range(len(rows)),rows,fontsize=7);ax.set_ylabel('Leiden cluster · r0.5')
   ax.set_xticks(range(len(page)),[ids_to_symbol[k] if ids_to_symbol[k] not in dup else ids_to_symbol[k]+' '+k for k in page],rotation=65,ha='right',fontsize=7)
   ax.set_title(p+' · '+name);save(fig,f'{name}_g{j//22+1}_c{h//36+1}')
matrixplot(top,'positive_dotplot');matrixplot(top,'positive_heatmap',True)
for k,ls in by_marker.items():matrixplot(ls,k+'_markers')
# Count object-level co-detection of multiple measured genes per family. This is an observation for manual review.
x,gg,obs=combine(p);label=z.set_index('cell_key').loc[obs.cell_key,'leiden_r0_5'].values
def present(fam):
 jj=np.flatnonzero(gg.feature_id.isin(by_marker[fam]));return x[:,jj].getnnz(1)>=2
epi=present('Epithelial');t=present('T_NK');my=present('Myeloid')
signals=pd.DataFrame({'cluster':label,'two_epi_and_two_immune':epi&(t|my),'two_epi_and_two_T_NK':epi&t,'two_epi_and_two_myeloid':epi&my})
review=template[['cluster','n_objects','dominant_sample','dominant_sample_fraction','median_nFeature','median_GEX_UMI','GEX_secondary_n']].merge(signals.groupby('cluster').mean(),on='cluster')
review['sample_bias_flag']=review.dominant_sample_fraction>=.85;review['low_feature_flag']=review.median_nFeature<np.quantile(z.nFeature,.25);review['co_detection_review_flag']=review.two_epi_and_two_immune>=.10
baseline=z['sample'].value_counts(normalize=True)
review['dominant_sample_platform_fraction']=review.dominant_sample.map(baseline)
review['dominant_sample_fraction_excess']=review.dominant_sample_fraction-review.dominant_sample_platform_fraction
write(review,d/'review_priorities.tsv')

