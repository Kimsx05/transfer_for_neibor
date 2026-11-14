from common import *
import importlib.util
sp=importlib.util.spec_from_file_location('plot',R/'code/13_replot.py');p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p);plt=p.plt;Normalize=p.Normalize
idx=pd.read_csv(R/'tables/figure_index.tsv',sep='\t').to_dict('records');p.index=idx
D=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');roi=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');bins=pd.read_parquet(R/'tables/Cervilla_bin8_measurement.parquet')
# All panels use direct saved raw or fitted outputs, not BANKSY features.
for r in roi[roi.roi.isin(['T_lineage_focus','vascular_focus','far_UMI_hotspot'])].itertuples():
 fig,axs=plt.subplots(1,4,figsize=(15,3.6),layout='constrained');specs=[];cb=[]
 for ax,key,label in zip(axs[:3],['Epi','Treg','Endo'],['Epi panel','Treg panel (not identity)','Endothelial panel']):
  q=bins[bins.x_um.between(r.x0,r.x1)&bins.y_um.between(r.y0,r.y1)].copy();q['raw_signal']=np.log1p(q[key+'_UMI']);vv,ss=p.spatial(ax,q,'raw_signal',label,lim=(r.x0,r.x1,r.y0,r.y1),vmax=4);bar=fig.colorbar(vv,ax=ax,shrink=.65);specs.append(ss);cb.append(bar.ax);ax.set_xlabel(label);ax.set_ylabel('')
 q=reg[reg.x_um.between(r.x0,r.x1)&reg.y_um.between(r.y0,r.y1)].copy();q['score']=q.MP10.where(q.score_eligible);vv,ss=p.spatial(axs[3],q,'score','Epi-only regional MP10',lim=(r.x0,r.x1,r.y0,r.y1),vmax=.14);axs[3].collections[-1].set_sizes([25]);bar=fig.colorbar(vv,ax=axs[3],shrink=.65);cb.append(bar.ax);specs.append(ss);axs[3].set_xlabel('Epi-only regional MP10');axs[3].set_ylabel('')
 for ax in axs:ax.set_xticks([]);ax.set_yticks([])
 p.save(fig,'ROI_targets_MP10_'+r.roi,specs,cb,'Raw log1p summed panel UMI on fixed8um bins; Treg signal is not validated identity; MP10 only eligible Epi-only region, gray unscorable')
# Target group detection dot panel (marker reuse is explicitly not independent truth).
e=pd.read_csv(R/'tables/target_group_marker_evidence.tsv',sep='\t');genes=['CD3D','CD3E','TRAC','FOXP3','IL2RA','CTLA4','IKZF2','CXCL13','TIGIT','PECAM1','VWF','CDH5','EMCN','EPCAM','KRT19','COL1A1'];groups=e.group.unique();z=e[e.gene.isin(genes)].copy();x=pd.Index(genes).get_indexer(z.gene);y=pd.Index(groups).get_indexer(z.group);fig,ax=plt.subplots(figsize=(11,4.5),layout='constrained');norm=Normalize(0,max(.2,z.mean_LogNormalize.quantile(.98)));sc=ax.scatter(x,y,s=300*z.positive_fraction,c=z.mean_LogNormalize,norm=norm,cmap='viridis');ax.set_xticks(range(len(genes)),genes,rotation=60,ha='right');ax.set_yticks(range(len(groups)),groups,fontsize=9);bar=fig.colorbar(sc,ax=ax,label='Mean LogNormalize expression');leg=ax.legend(handles=[ax.scatter([],[],s=300*t,color='gray',label=f'{t:.0%}')for t in [.1,.5,1]],loc='upper left',bbox_to_anchor=(1.25,1),frameon=False,title='Detected');p.save(fig,'target_group_marker_support',[('continuous',norm,'viridis','Mean LogNormalize'),('categorical',{'10% detected':'#bbbbbb','50% detected':'#888888','100% detected':'#555555'})],[bar.ax,leg],'Dot area=detected fraction; candidate Treg group is not a validated identity; TIGIT excluded from reference fitting')
# Fixed-ROI fitted quantity comparison: units have different definitions, no interpolation.
units={}
for u in ['bin8','bin16']:
 q=pd.read_parquet(R/'tables'/f'Cervilla_{u}_measured.parquet');w=pd.read_parquet(R/'tables'/f'{u}_ROI_contributions.parquet');q=q.merge(w,on='bin_id',validate='one_to_one');q['Epi_support']=q.contribution_Epithelial.where(q.reliable_fit);ec=[c for c in q if c.startswith('contribution_')and ('Endothelial' in c or c.endswith(' EC'))];q['Endo_support']=q[ec].sum(1).where(q.reliable_fit);q['MP10_show']=q.mp10_historical_auc.where(q.reliable_fit&(q.Epi_support>=.8)&~((q.Epi_detected>=2)&((q.Fib_detected>=2)|(q.T_detected>=2)|(q.Endo_detected>=2)))&(q.MP10_detected>=5)&(q.nFeature_RNA>=905));units[u]=q
units['cell']=D.assign(Epi_support=D.type_supported.eq('Epithelial').astype(float),Endo_support=D.type_supported.eq('Endothelial').astype(float),MP10_show=D.mp10_historical_auc.where(D.mp10_single_score_eligible))
for r in roi[roi.roi.isin(['T_lineage_focus','vascular_focus','far_UMI_hotspot'])].itertuples():
 fig,axs=plt.subplots(3,3,figsize=(10,10),layout='constrained');specs=[];cb=[]
 for i,(u,q)in enumerate(units.items()):
  z=q[q.x_um.between(r.x0,r.x1)&q.y_um.between(r.y0,r.y1)]
  for j,(key,lab,mx)in enumerate([('Epi_support','Epi contribution / cell indicator',1),('Endo_support','Endo contribution / cell indicator',1),('MP10_show','Single-unit MP10 (eligible only)',.14)]):
   vv,ss=p.spatial(axs[i,j],z,key,lab,lim=(r.x0,r.x1,r.y0,r.y1),vmax=mx);bar=fig.colorbar(vv,ax=axs[i,j],shrink=.65);cb.append(bar.ax);specs.append(ss);axs[i,j].set_xlabel(lab,fontsize=9);axs[i,j].set_ylabel(u);axs[i,j].set_xticks([]);axs[i,j].set_yticks([])
 p.save(fig,'ROI_model_and_score_'+r.roi,specs[:3],cb,'Rows bin8/bin16/cell; model contribution versus binary supported cell identity are distinct. Gray=not reliable/scorable; no fraction division; regional score is main, single-unit scores diagnostic')
# D1 light registration review at two inherited ROIs. No fitting or spatial relationship tests.
assets=json.loads((A/'segmentation_assignment_v2/D1_Johnson/metadata/assets.json').read_text());tr=json.loads((A/'segmentation_assignment_v2/D1_Johnson/metadata/transform.json').read_text());im=p.Image.open(H/assets['image_project_relative']);mpp=tr['native_mpp'];ox,oy=tr['fullres_crop_xyxy'][:2];dd=pd.read_parquet(R/'tables/D1_Johnson_bin8_measurement.parquet');dr=[]
for name,x0,y0 in [('conflict_focus',7039.901,11967.634),('reference_agree_control',8704.,13312.)]:
 x1=x0+256;y1=y0+256;dr.append(dict(section='D1_Johnson',roi=name,x0=x0,x1=x1,y0=y0,y1=y1,reason='inherited ROI, light target/registration check only'));fig,axs=plt.subplots(1,4,figsize=(14,3.7),layout='constrained');crop=im.crop((x0/mpp-ox,y0/mpp-oy,x1/mpp-ox,y1/mpp-oy));axs[0].imshow(crop,extent=[x0,x1,y1,y0]);axs[0].set_xlabel('H&E (crop origin preserved)');specs=[];cb=[]
 for ax,key in zip(axs[1:],['Epi','T','Endo']):
  z=dd[dd.x_um.between(x0,x1)&dd.y_um.between(y0,y1)].copy();z['signal']=np.log1p(z[key+'_UMI']);vv,ss=p.spatial(ax,z,'signal','log1p raw '+key+' panel',lim=(x0,x1,y0,y1),vmax=4);bar=fig.colorbar(vv,ax=ax,shrink=.65);cb.append(bar.ax);specs.append(ss);ax.set_xlabel(key+' raw panel');ax.set_ylabel('')
 for ax in axs:ax.set_xlim(x0,x1);ax.set_ylim(y1,y0);ax.set_xticks([]);ax.set_yticks([])
 axs[0].plot([x0+10,x0+60],[y1-10]*2,'k-',lw=2);axs[0].text(x0+35,y1-15,'50 μm',ha='center');p.save(fig,'D1_light_'+name,specs,cb,'Existing D1 crop and physical coordinates preserved; raw panel evidence, not new fine labels or associations')
put(pd.DataFrame(dr),'D1_light_roi_manifest.tsv')
# Figure genuinely generated from the R loader readback, not an earlier cached input.
f=R/'tables/R_loader_readback.tsv.gz'
if f.exists():
 q=pd.read_csv(f,sep='\t');p.catplot(q,'r50_niche','R_loader_readback_niche50',note='Actual reread original Seurat counts plus ID-aligned new metadata; no new full expression object saved')
# Compact stability and historical measurement diagnostics.
a=pd.read_csv(R/'tables/niche_candidate_metrics.tsv',sep='\t');fig,axs=plt.subplots(1,2,figsize=(9,4),layout='constrained');specs=[];legs=[]
for ax,r in zip(axs,[50,200]):
 q=a[a.radius_um==r]
 for key,col in [('seed_ARI_min','#1b9e77'),('block_ARI','#d95f02'),('block_min_cluster_Jaccard','#7570b3')]:ax.plot(q.K,q[key],'-o',c=col,label=key)
 ax.set_xlabel(f'K ({r} μm)');ax.set_ylabel('Agreement (not accuracy)');ax.set_ylim(0,1.03);legs.append(ax.legend(fontsize=8,frameon=False))
p.save(fig,'niche_stability',[('categorical',{'Seed ARI':'#1b9e77','Block ARI':'#d95f02','Min matched Jaccard':'#7570b3'})],legs,'Fixed seeds and one spatial-block refit; no MP10 selection')
fig,axs=plt.subplots(1,2,figsize=(9,4.5),layout='constrained');q=reg[reg.score_eligible]
for ax,col in zip(axs,['MP10_half','MP10_depth3000']):ax.scatter(q.MP10,q[col],s=6,color='#555555',alpha=.4,rasterized=True);ax.set_xlabel('Historical MP10, original depth');ax.set_ylabel(col)
p.save(fig,'MP10_historical_depth_sensitivity',note='Fixed binomial thinning and historical mean-tie ranks; no target-association based score selection')
put(pd.DataFrame(p.index),'figure_index.tsv');(R/'FIGURE_INDEX.md').write_text('# Figure index\n\n'+'\n'.join('- **'+x['figure']+'**: '+x['note'] for x in p.index)+'\n\nAll final figures: 300 dpi PNG and PDF. H&E and dense points raster; labels/axes/legends vector. Non-evaluable gray is explicitly defined. Preview PNGs are raw measurement checks.\n')
print('SUPPLEMENT_PLOTS_COMPLETE')
