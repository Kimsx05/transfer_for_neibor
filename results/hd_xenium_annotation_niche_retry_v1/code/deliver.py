from base_data import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from scipy.spatial import cKDTree
from collections import Counter
import colorsys,hashlib
pal=json.loads((V/'shared/colors.json').read_text());pal.update(Mixed='#cd853f',major_unresolved='#a9a9a9',major_only='#4477aa',mixed_or_suspect='#cc8844',fine_supported='#228833',unassigned_QC='#dadada')
for k in range(8):
 for platform in ['hd','xenium']:pal[f'{platform}_N{k+1}']=plt.matplotlib.colors.to_hex(plt.cm.Set2(k/7))
def color(k):
 if k not in pal:
  h=int(hashlib.sha256(k.encode()).hexdigest()[:6],16);pal[k]=plt.matplotlib.colors.to_hex(colorsys.hls_to_rgb((h%360)/360,.49,.61))
 return pal[k]
def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(O/'figures'/(name+'.'+ext),bbox_inches='tight',dpi=155)
 plt.close(fig)
def spatial(z,col,title,name,roi=False):
 fig,ax=plt.subplots(figsize=(6.2,4.8));leg=[];values=z[col].astype(str)
 for k in sorted(values.unique()):
  sub=z[values==k];cc=color(k);ax.scatter(sub.x_um,sub.y_um,s=4 if roi else .34,c=cc,linewidth=0,rasterized=True);leg.append(Line2D([],[],color=cc,marker='o',ls='',ms=4,label=k))
 ax.legend(handles=leg,loc='center left',bbox_to_anchor=(1, .5),frameon=False,fontsize=6,ncol=2 if len(leg)>14 else 1);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_title(title,fontsize=10);save(fig,name)
platform=sys.argv[1];samples=XE if platform=='xenium' else HD
allobs=[];dots=[];techn=[];cover=[];roirows=[];geom=[];corr=[]
genes=list(dict.fromkeys(sum(M.values(),[])+sum((str(v).split() for v in tax.anchors),[])))
for s in samples:
 d=paths(s);z=pd.read_parquet(d/'annotation.parquet').join(pd.read_parquet(d/'niche.parquet'),validate='one_to_one');a=dataset(s);assert z.index.equals(a.obs_names) and z.niche_id.notna().all();z.index.name='original_cell_id';z.to_parquet(d/'analysis_index.parquet');allobs.append(z)
 for col in ['major_identity','display_label','review_level','niche_id','conflict']:
  spatial(z,col,s+' · '+col,s+'_'+col)
 if platform=='hd':spatial(z,'source',s+' · reconstruction source',s+'_source')
 # Representative locations selected by group density with a minimum denominator, independent of target states.
 tiles=pd.DataFrame({'tx':np.floor(z.x_um/140),'ty':np.floor(z.y_um/140)},index=z.index)
 for name,signal in [('supported',z.review_level.eq('fine_supported')),('conflict',z.review_level.eq('mixed_or_suspect')),('unresolved',z.review_level.eq('major_unresolved'))]+([('secondary',z.source.eq('GEX_secondary'))] if platform=='hd' else []):
  tab=tiles.assign(target=signal).groupby(['tx','ty']).target.agg(['sum','size']);tab=tab[tab['size']>=15];k=tab['sum'].idxmax();idx=(tiles.tx==k[0])&(tiles.ty==k[1]);sub=z[idx]
  for col in ['niche_id','display_label']+(['source'] if name=='secondary' else []):spatial(sub,col,s+' · '+name+' ROI',s+'_ROI_'+name+'_'+col,True)
  roirows.append((s,name,*[int(v*140) for v in k],len(sub),int(signal[idx].sum())))
  # Raw detected marker maps in the same small window.
  marker={'supported':'EPCAM','conflict':'PTPRC','unresolved':'KRT8','secondary':'KRT8'}[name]
  if marker in a.var_names:
   raw=a.X[idx.to_numpy(),a.var_names.get_loc(marker)].toarray().ravel();fig,ax=plt.subplots(figsize=(4,4));p=ax.scatter(sub.x_um,sub.y_um,c=np.log1p(raw),s=9,linewidth=0,cmap='viridis');fig.colorbar(p,ax=ax,label='log(1 + raw counts)');ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_title(s+' · '+name+' · '+marker);save(fig,s+'_ROI_'+name+'_raw')
 g=pd.Index(genes).intersection(a.var_names,sort=False);xx=a.X[:,a.var_names.get_indexer(g)];library=np.maximum(1,a.X.sum(1).A1);yy=xx.astype('f4').multiply((1e4/library)[:,None]).tocsr();yy.data=np.log1p(yy.data)
 for k,ids in z.groupby('lineage_cluster').indices.items():
  n=len(ids);cnt=xx[ids].getnnz(0);ss=yy[ids].sum(0).A1;dots.extend(zip([s]*len(g),[k]*len(g),g,[n]*len(g),cnt,ss))
 for niche,unit in z.groupby('niche_id'):
  techn.append((s,niche,len(unit),len(unit)/len(z),unit.n_counts.median(),unit.nFeature.median(),unit.area_um2.median(),unit.source.eq('GEX_secondary').mean(),unit.geometry_flag.mean(),unit.review_level.eq('mixed_or_suspect').mean(),unit.review_level.eq('major_unresolved').mean()))
 lat=np.load(d/'latent.npy');idx=np.isfinite(lat[:,0]);q=z[idx]
 for feature in ['n_counts','nFeature','area_um2']:
  for j in range(lat.shape[1]):corr.append((s,feature,j,float(spearmanr(lat[idx,j],q[feature]).statistic)))
 for k,sub in z.groupby('source'):
  geom.append((s,k,len(sub),sub.n_counts.median(),sub.nFeature.median(),sub.area_um2.median(),sub.geometry_flag.mean(),sub.review_level.eq('mixed_or_suspect').mean(),sub.major_identity.eq('Unknown').mean(),sub.eligible_expression.mean()))
 if platform=='hd':
  rtab=pd.read_csv(O/'tables'/(s+'_ROI.tsv'),sep='\t');members=pd.read_parquet(d/'members.parquet',columns=['array_row','array_col','final_label','cell_label']);idmap=pd.Series(z.index,index=a.obs.label)
  for r in rtab.itertuples():
   x=members[members.array_row.between(r.row,r.row+83)&members.array_col.between(r.col,r.col+83)]
   lab=x.final_label if r.route=='joint' else x.cell_label;sub=z.loc[idmap.loc[np.unique(lab[lab>0])]]
   cover.append((s,r.roi,r.route,len(sub),sub.n_counts.median(),sub.nFeature.median(),sub.area_um2.median(),sub.review_level.eq('mixed_or_suspect').mean(),int(sub.geometry_flag.sum()),r.unassigned_umi))
 del a
allz=pd.concat(allobs);write(pd.DataFrame(corr,columns=['sample','technical','latent_dimension','rho']),platform+'_embedding_QC.tsv')
write(pd.DataFrame(roirows,columns=['sample','rationale','x0_um','y0_um','n','selected_property_n']),platform+'_review_ROI.tsv')
write(pd.DataFrame(geom,columns=['sample','source','n','median_UMI','median_genes','median_area','geometry_flag_fraction','mixed_fraction','unknown_fraction','expression_eligible_fraction']),platform+'_source_QC.tsv')
if cover:write(pd.DataFrame(cover,columns=['sample','roi','route','n_touching_core','median_umi','median_genes','median_area','mixed_fraction','suspicious_geometry','unassigned_core_umi']),'hd_ROI_comparison.tsv')
write(pd.DataFrame(techn,columns=['sample','niche','n','sample_fraction','median_counts','median_genes','median_area','secondary_fraction','geometry_flag_fraction','mixed_fraction','unknown_fraction']),platform+'_niche_support.tsv')
ct=allz.groupby(['sample','niche_id','major_identity','display_label']).size().rename('n').reset_index();ct['within_niche_fraction']=ct.n/ct.groupby(['sample','niche_id']).n.transform('sum');write(ct,platform+'_niche_composition.tsv')
dd=pd.DataFrame(dots,columns=['sample','cluster','gene','n','detected','sum_log_cp10k']);dd.to_parquet(O/'tables'/(platform+'_cluster_markers.parquet'),index=False)
for group in ['T/NK','Myeloid','Endothelial','Fibroblast','Epithelial']:
 wanted=list(dict.fromkeys(M.get(group,[])+sum([str(v).split()[:3] for v in tax[tax.major==group].anchors],[])))[:28]
 sub=dd[dd.cluster.str.startswith(group+'.')&dd.gene.isin(wanted)].groupby(['cluster','gene'])[['n','detected','sum_log_cp10k']].sum().reset_index()
 if sub.empty:continue
 fig,ax=plt.subplots(figsize=(max(6,len(wanted)*.32),max(3,sub.cluster.nunique()*.26)));cats=sorted(sub.cluster.unique());gx=[g for g in wanted if g in set(sub.gene)]
 pts=ax.scatter(sub.gene.map(dict(zip(gx,range(len(gx))))),sub.cluster.map(dict(zip(cats,range(len(cats))))),s=1+100*sub.detected/sub.n,c=sub.sum_log_cp10k/sub.n,cmap='viridis');ax.set_xticks(range(len(gx)),gx,rotation=65,ha='right',fontsize=7);ax.set_yticks(range(len(cats)),cats,fontsize=7);fig.colorbar(pts,ax=ax,label='mean log-normalized counts');ax.set_title(platform+' · '+group+' · size = fraction detected');save(fig,platform+'_dot_'+group.replace('/','_'))
(O/'shared'/('colors_'+platform+'.json')).write_text(json.dumps(pal,indent=2))

