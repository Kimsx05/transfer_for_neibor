"""Render only saved measurement/annotation/model-result tables and linked H&E/masks. No inference."""
from common import *
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt,matplotlib.font_manager as fm
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
fm.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1,'lines.linewidth':1.2,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
COL={'0':'#bdbdbd','Epithelial':'#f5b264','Fibroblast':'#8c6bb1','Mural':'#756bb1','Endothelial':'#d95f02','Myeloid':'#80a2d3','T_subtype_unresolved':'#57c1b5','Treg':'#238b45','B':'#bbcf5c','Plasma':'#82cfad','Mast':'#00baee','NK':'#1b9e77','Unknown':'#bdbdbd'}
index=[]
def save(fig,name,specs=None,axes_to_hide=None,note=''):
 specs=specs or [];axes_to_hide=axes_to_hide or [];fig.savefig(R/'figures'/f'{name}.png',dpi=300);fig.savefig(R/'figures'/f'{name}.pdf')
 if specs:
  for a in axes_to_hide:a.set_visible(False)
  fig.savefig(R/'figures'/f'{name}_panel_only.pdf')
  for a in axes_to_hide:a.set_visible(True)
  lf,aa=plt.subplots(len(specs),1,figsize=(8,max(2,1.8*len(specs))),squeeze=False,layout='constrained')
  for ax,sp in zip(aa[:,0],specs):
   if sp[0]=='continuous':lf.colorbar(plt.cm.ScalarMappable(norm=sp[1],cmap=sp[2]),cax=ax,orientation='horizontal',label=sp[3])
   else:ax.axis('off');ax.legend(handles=[Patch(color=v,label=k)for k,v in sp[1].items()],loc='center',frameon=False,ncol=1,fontsize=10)
  lf.savefig(R/'figures'/f'{name}_legend_only.pdf');plt.close(lf)
 plt.close(fig);index.append(dict(figure=name,note=note,source='tables/*.parquet and TSV; recorded ROI/linked images',format='white Arial12; 300dpi PNG; vector PDF text/legends with raster dense points/H&E'))
def spatial(ax,d,col=None,label=None,lim=None,cmap='viridis',vmin=None,vmax=None):
 if col is None:ax.scatter(d.x_um,d.y_um,s=.5,color='#ededed',rasterized=True);return None
 vals=d[col].to_numpy();ok=np.isfinite(vals);vmax=vmax if vmax is not None else (np.quantile(vals[ok],.99) if ok.any() else 1);vmin=0 if vmin is None else vmin;vmax=max(vmax,vmin+1e-9);norm=Normalize(vmin,vmax);ax.scatter(d.x_um,d.y_um,s=1,color='#d9d9d9',rasterized=True);z=ax.scatter(d.x_um[ok],d.y_um[ok],c=vals[ok],s=3,cmap=cmap,norm=norm,rasterized=True);ax.set_aspect('equal');ax.set_xlabel('x (μm)');ax.set_ylabel('y (μm)');ax.invert_yaxis()
 if lim:ax.set_xlim(lim[:2]);ax.set_ylim(lim[3],lim[2])
 return z,('continuous',norm,cmap,label or col)
def catplot(d,col,name,lim=None,note=''):
 cats=sorted(d[col].dropna().unique());palette={str(k):COL.get(str(k),matplotlib.colors.to_hex(plt.cm.tab20(i%20)))for i,k in enumerate(cats)};fig,ax=plt.subplots(figsize=(7.5,5.5),layout='constrained')
 for k in cats:
  q=d[d[col]==k];ax.scatter(q.x_um,q.y_um,c=palette[str(k)],s=1,rasterized=True)
 ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (μm)');ax.set_ylabel('y (μm)');leg=ax.legend(handles=[Patch(color=v,label=k)for k,v in palette.items()],bbox_to_anchor=(1,1),loc='upper left',frameon=False,fontsize=9)
 if lim:ax.set_xlim(lim[:2]);ax.set_ylim(lim[3],lim[2])
 save(fig,name,[('categorical',palette)],[leg],note)

def main():
 d=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');bins=pd.read_parquet(R/'tables/Cervilla_bin8_measurement.parquet');lim=(bins.x_um.min()-80,bins.x_um.max()+80,bins.y_um.min()-80,bins.y_um.max()+80)
 catplot(d,'type_supported','cell_supported_lineages',lim,'Candidate lineage support, not validated cell identities; Unknown retained; missing supplied image area has no data')
 for name,q,key,lab in [('MP10_epithelial_regions',reg.assign(MP10_plot=reg.MP10.where(reg.score_eligible)),'MP10_plot','Epi-only MP10 rank-AUC'),('cell_detection',d,'nFeature_RNA','Detected genes'),('cell_Endothelial_raw',d,'Endo_UMI','Endothelial-panel raw UMI')]:
  fig,ax=plt.subplots(figsize=(6.8,5.4),layout='constrained');v,sp=spatial(ax,q,key,lab,lim);cb=fig.colorbar(v,ax=ax,label=lab);save(fig,name,[sp],[cb.ax],'Original RNA evidence; gray = not scorable for MP10, white outside supplied support is unmeasured')
 # Explicit Treg evaluability figure, not an abundance map.
 dd=d.copy();dd['Treg_diagnostic']=np.where(dd.fine_prediction.eq('Treg'),'Unvalidated Treg candidate','Other QC objects');catplot(dd,'Treg_diagnostic','Treg_candidate_not_confirmed',lim,'No reliable Treg population established; candidates cannot be interpreted as Treg abundance')
 # Registration overview, display image under supplied-bin footprint and reference ROIs.
 assets=json.loads((A/'segmentation_assignment_v2/Cervilla/metadata/assets.json').read_text());im=Image.open(H/assets['image_project_relative']);im.thumbnail((2200,2200));ro=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');fig,ax=plt.subplots(figsize=(7.2,5.7),layout='constrained');ax.imshow(im,extent=[0,11418*.8636122748145504,8964*.8636122748145504,0]);ax.scatter(bins.x_um,bins.y_um,s=.05,c='#aaaaaa',alpha=.07,rasterized=True)
 for r in ro.itertuples():ax.add_patch(plt.Rectangle((r.x0,r.y0),r.x1-r.x0,r.y1-r.y0,fill=False,edgecolor='cyan',lw=1.1));ax.text(r.x0,r.y0,r.roi,fontsize=8,color='black')
 ax.set_xlabel('x (μm)');ax.set_ylabel('y (μm)');save(fig,'HE_ROI_registration',note='Loupe-recovered H&E embedded as raster; image outside supplied footprint is unmeasured, not negative')
 # Bin/cell raw evidence and boundary overlay; no RNA smoothing.
 imagefull=Image.open(H/assets['image_project_relative']);maskn=sparse.load_npz(A/'segmentation_assignment_v2/Cervilla/inputs/nuclei.npz').tocsr();maskc=sparse.load_npz(A/'segmentation_assignment_v2/Cervilla/inputs/assignment_mask.npz').tocsr();tr=json.loads((A/'segmentation_assignment_v2/Cervilla/metadata/transform.json').read_text());mpp=tr['native_mpp'];sx,sy=tr['local_to_model_scale_xy'];b16=pd.read_parquet(R/'tables/Cervilla_bin16_measurement.parquet')
 from skimage.segmentation import find_boundaries
 for r in ro.itertuples():
  fig,axs=plt.subplots(1,4,figsize=(15,3.5),layout='constrained');specs=[];cbs=[];crop=imagefull.crop((int(r.x0/mpp),int(r.y0/mpp),int(r.x1/mpp),int(r.y1/mpp)));axs[0].imshow(crop,extent=[r.x0,r.x1,r.y1,r.y0]);xx0=int((r.x0/mpp+.5)*sx-.5);xx1=int((r.x1/mpp+.5)*sx-.5);yy0=int((r.y0/mpp+.5)*sy-.5);yy1=int((r.y1/mpp+.5)*sy-.5)
  for mask,color in [(maskc,'#ffcc00'),(maskn,'#00ffff')]:
   sub=mask[yy0:yy1,xx0:xx1].toarray();bound=find_boundaries(sub,mode='inner');rgba=np.zeros((*bound.shape,4));rgba[bound]=matplotlib.colors.to_rgba(color,alpha=.7);axs[0].imshow(rgba,extent=[((xx0+.5)/sx-.5)*mpp,((xx1+.5)/sx-.5)*mpp,((yy1+.5)/sy-.5)*mpp,((yy0+.5)/sy-.5)*mpp])
  axs[0].set_xlabel('H&E / nucleus / 8μm extent')
  for ax,q,key,title in zip(axs[1:],[bins,b16,d],['T_UMI','T_UMI','T_UMI'],['8μm bin raw T panel','16μm bin raw T panel','Cell raw T panel']):
   z=q[q.x_um.between(r.x0,r.x1)&q.y_um.between(r.y0,r.y1)].copy();z['signal']=np.log1p(z[key]);vv,sp=spatial(ax,z,'signal','log1p raw T-panel UMI',vmax=3);cb=fig.colorbar(vv,ax=ax,shrink=.7);specs.append(sp);cbs.append(cb.ax);ax.set_xlabel(title);ax.set_ylabel('')
  for ax in axs:ax.set_xlim(r.x0,r.x1);ax.set_ylim(r.y1,r.y0);ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal')
  axs[0].plot([r.x0+10,r.x0+60],[r.y1-10]*2,'k-',lw=2);axs[0].text(r.x0+35,r.y1-15,'50 μm',ha='center',fontsize=9);save(fig,'ROI_bin_cell_'+r.roi,specs,cbs,'Same physical ROI; nucleus cyan and inferred extent yellow; raw T panel identical scale; original bin16 grid includes all supplied counts')
 del maskn,maskc
 # Score measurement robustness and depth; direct saved data only.
 for unit in ['cell','bin8','bin16']:
  v=pd.read_csv(R/'tables'/f'{unit}_score_seed_validation.tsv.gz',sep='\t');fig,ax=plt.subplots(figsize=(5,4.7),layout='constrained')
  for mode,col in [('auc_original','#d95f02'),('auc_observed','#1b9e77')]:
   z=v.pivot(index='id',columns='seed',values=mode);ax.scatter(z[123],z[124],s=9,c=col,label={'auc_original':'Native random zero ties','auc_observed':'Zero ranks excluded'}[mode],alpha=.5)
  ax.set_xlabel('Seed 123 AUCell');ax.set_ylabel('Seed 124 AUCell');leg=ax.legend(frameon=False);save(fig,'score_ties_'+unit,[('categorical',{'Native AUCell random zeros':'#d95f02','zero ranks excluded':'#1b9e77'})],[leg],'Fixed detection-stratified measurement subset; not an association test')
 for r in [50,200]:
  n=pd.read_parquet(R/'tables'/f'Cervilla_cell_neighborhood_r{r}.parquet');catplot(n,'niche',f'niche_r{r}',lim,'Labels only on Epi centers; 0 is unresolved; full known composition includes Epi')
  for key,label in [('known_fraction','Known objects / all QC neighbors'),('observed_area_fraction','Observed area / circle area')]:
   fig,ax=plt.subplots(figsize=(6.8,5.4),layout='constrained');v,sp=spatial(ax,n,key,label,lim,vmax=1);cb=fig.colorbar(v,ax=ax,label=label);save(fig,f'niche_r{r}_{key}',[sp],[cb.ax],'Quality field never used as clustering feature')
  comp=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t');comp=comp[comp.radius_um==r];types=[c[2:]for c in n if c.startswith('P_')];fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained');base=np.zeros(len(comp))
  palette={t:COL.get(t,'#969696')for t in types}
  for t in types:ax.bar(comp.niche.astype(str),comp[t],bottom=base,color=palette[t],label=t,width=.8);base+=comp[t]
  ax.set_ylabel('Known-type relative composition');ax.set_xlabel(f'Niche ({r} μm)');leg=ax.legend(bbox_to_anchor=(1,1),loc='upper left',frameon=False,fontsize=9);save(fig,f'niche_composition_r{r}',[('categorical',palette)],[leg],'Neighbor object-count proportions; Unknown excluded from conditional denominator')
  v=pd.read_parquet(R/'tables'/f'Cervilla_region_association_r{r}.parquet');v=v[v.association_eligible]
  for key,label in [('Endothelial_density','Supported endothelial objects / observed mm²'),('Endothelial_distance','Distance to endothelial support (μm)')]:
   q=v.dropna(subset=[key,'MP10']);fig,ax=plt.subplots(figsize=(5.5,4.7),layout='constrained');ax.scatter(q.MP10,q[key],s=7,color='#505050',alpha=.4,rasterized=True);ax.set_xlabel('Epi-only MP10 rank-AUC');ax.set_ylabel(label);save(fig,f'MP10_{key}_r{r}',note='One point per fixed Epi-only region; no bin-level pseudoreplication; censoring retained in table')
  q=v[v.niche_region_summary>0];groups=sorted(q.niche_region_summary.unique());fig,ax=plt.subplots(figsize=(6.5,4.6),layout='constrained');ax.boxplot([q.loc[q.niche_region_summary==k,'MP10']for k in groups],tick_labels=[str(k)for k in groups],showfliers=False);ax.set_xlabel(f'Majority pre-frozen niche ({r} μm)');ax.set_ylabel('Epi-only MP10 rank-AUC');save(fig,f'niche_MP10_r{r}',note='Region summary requires >=60% of its frozen Epi centers in one niche; mixed regions remain explicit')
 a=pd.read_csv(R/'tables/association_results.tsv',sep='\t');q=a[a.P.notna()&a.family.eq('continuous_primary')];fig,ax=plt.subplots(figsize=(7.2,5),layout='constrained');yy=np.arange(len(q));ax.errorbar(q.effect,yy,xerr=np.vstack([q.effect-q.low,q.high-q.effect]),fmt='o',color='black',capsize=3);ax.axvline(0,color='#999999',lw=1);ax.set_yticks(yy,[f'{r.radius_um} μm {r.predictor} / {r.model}'for r in q.itertuples()],fontsize=9);ax.set_xlabel('Standardized rank effect (spatial-block 95% CI)');save(fig,'association_effect_intervals',note='499 block bootstrap CI; BH q and restricted-null block-sign P in table; single section')
 put(pd.DataFrame(index),'figure_index.tsv');(R/'FIGURE_INDEX.md').write_text('# Figure index\n\n'+'\n'.join('- **'+x['figure']+'**: '+x['note'] for x in index)+'\n\nH&E and dense points are raster; labels, coordinates, and legends are vector. Gray status definitions differ by plot and are explicitly recorded. Preview figures are measurements only; final figures provide PDF/PNG and separate legends.\n')
if __name__=='__main__':main()
