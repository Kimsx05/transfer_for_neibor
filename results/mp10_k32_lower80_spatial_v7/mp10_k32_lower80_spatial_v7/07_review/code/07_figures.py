from common import *
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from scipy.stats import spearmanr
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
blue='#173B75';cblue=LinearSegmentedColormap.from_list('k32',['#f0f3f8',blue]);grey=LinearSegmentedColormap.from_list('other',['#fafafa','#555555']);F=R/'06_figures';F.mkdir(exist_ok=True);captions=[]
def save(fig,p,data,cap):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(p)+'.pdf',bbox_inches='tight',transparent=True);fig.savefig(str(p)+'.png',dpi=300,bbox_inches='tight',transparent=True);plt.close(fig);wt(data,str(p)+'_plot_data.tsv.gz');Path(str(p)+'_caption.txt').write_text(cap+'\n');captions.append({'figure':str(p.relative_to(R)),'caption':cap})
def scattermap(ax,d,col,cmap='viridis',vmax=None,zero_black=False):
 x=d[col].to_numpy(float);ok=np.isfinite(x);ax.scatter(d.x,d.y,c='#bdbdbd',s=8,linewidths=0);positive=ok&(x!=0) if zero_black else ok
 if zero_black:ax.scatter(d.x[ok&(x==0)],d.y[ok&(x==0)],c='black',s=8,linewidths=0)
 im=ax.scatter(d.x[positive],d.y[positive],c=x[positive],s=8,cmap=cmap,norm=Normalize(0,vmax if vmax is not None and vmax>0 else max(np.nanmax(x),1e-8)),linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('Measured spatial coordinates');ax.figure.colorbar(im,ax=ax,shrink=.65,ticks=np.linspace(0,vmax,5) if vmax is not None and vmax>0 else None,format='%.3g')
 from matplotlib.patches import Patch
 handles=([Patch(color='black',label='True score zero')] if zero_black else [])+([Patch(color='#bdbdbd',label='Missing / unstable')] if (~ok).any() else [])
 if handles:ax.legend(handles=handles,loc='upper left',fontsize=5,frameon=False)
 return im
# Common abundance display limits from available complete posterior tables; raw values retained.
spatial={};statuses={}
for sec in config()['sections']:
 o=outmap(sec);p=o/'raw_expression_plot_data.tsv.gz'
 if not p.exists():continue
 d=pd.read_csv(p,sep='\t');status='RAW_EXPRESSION_ONLY'
 if (o/'RESULT_STATUS.json').exists():
  st=json.loads((o/'RESULT_STATUS.json').read_text());ab=pd.read_csv(o/'abundance_mean.tsv.gz',sep='\t',index_col=0).reindex(d.spot_id);d['K32_mean']=ab[K].to_numpy();d['lower80_mean']=ab[L].to_numpy();d['vCAF_mean']=ab.vCAF.to_numpy();q=pd.read_csv(o/'abundance_q05.tsv.gz',sep='\t',index_col=0);d['K32_q05']=q.loc[d.spot_id,K].to_numpy();der=pd.read_csv(o/'posterior_derived.tsv.gz',sep='\t');d=d.merge(der,on='spot_id',validate='one_to_one');status=st['numerical_status']+'; epochs='+str(st['epochs'])
 spatial[sec]=d;statuses[sec]=status
limits={c:float(np.nanquantile(np.concatenate([d[c].to_numpy() for d in spatial.values() if c in d]),.99)) for c in ['K32_mean','lower80_mean','Epi_total_mean','vCAF_mean','K32_q05'] if any(c in d for d in spatial.values())};clip=[];spcorr=[]
for sec,d in spatial.items():
 o=F/'by_dataset'/sec.split('_')[0]/sec;display=pd.read_csv(R/'00_manifest/DISPLAY_GENES.tsv',sep='\t').gene.tolist();rawcols=['MP10_Top20_UCell','MP10_Top20_AUCell','EPCAM_log1pCP10k','KRT8_log1pCP10k','KRT18_log1pCP10k']+[g+'_log1pCP10k' for g in display]
 fig,axs=plt.subplots(2,4,figsize=(15,8 if np.ptp(d.y)>np.ptp(d.x) else 5.5));
 for ax,col in zip(axs.flat,rawcols):
  if col in d and d[col].notna().any():scattermap(ax,d,col,'plasma' if col.startswith('MP10') else 'viridis',zero_black=col.startswith('MP10'));ax.set_title(col.replace('_log1pCP10k',''))
  else:ax.text(.5,.5,col+'\nNot measured',ha='center',va='center',transform=ax.transAxes);ax.set_axis_off()
 fig.suptitle(sec+' | full measured RNA');save(fig,o/'raw_expression_overview',d[['spot_id','x','y','nCount_full']+rawcols],'Frozen MP10 Top20 scored in full measured RNA universe. Rank ties averaged; maxRank1500; AUCell cutoff ceil(0.05*n_genes). Exact zero black, missing/unmeasured grey. Gene expression log1p(count/total RNA*10000). Shared program expression is consistency, not independent validation. Coordinate-only view; no unverified histology overlay.')
 fig,ax=plt.subplots(figsize=(6,5));ax.scatter(d.x,d.y,s=8,c='#888888');ax.invert_yaxis();ax.set_aspect('equal');ax.set_xlabel('x');ax.set_ylabel('y');save(fig,o/'coordinate_basemap',d[['spot_id','x','y']],'Original QC-retained spot coordinates, same orientation and spot set as every panel; no fabricated histology.')
 if 'K32_mean' not in d:continue
 cols=['K32_mean','lower80_mean','Epi_total_mean','I_state_mean','K32_q05','vCAF_mean'];fig,axs=plt.subplots(2,3,figsize=(12,8 if np.ptp(d.y)>np.ptp(d.x) else 5.2))
 for ax,col in zip(axs.flat,cols):
  vm=1 if col=='I_state_mean' else limits[col];scattermap(ax,d,col,cblue if col.startswith('K32') or col=='I_state_mean' else grey if col=='lower80_mean' else 'viridis',vm);ax.set_title(col)
  clip.append({'section':sec,'variable':col,'vmax':vm,'n_above_display_max':int((d[col]>vm).sum()),'n_missing':int(d[col].isna().sum()),'scale':'0-1 fixed' if col=='I_state_mean' else 'pooled available sections q99'})
 fig.suptitle(sec+' | '+statuses[sec]);cap='Exploratory 31-type mapping; model-assigned abundance is not a validated physical cell count. '+statuses[sec]+'. Shared abundance scales use pooled section q99 (clipping counts in display_limits.tsv; full values retained); I_state uses draw-level ratio then posterior mean with fixed0–1. Grey ratio means >5% draws with Epi_total<0.1. q05 is a lower abundance quantile, not positivity probability. vCAF retained as competing type.';save(fig,o/'mapping_overview',d,cap)
 for col in cols:
  fig,ax=plt.subplots(figsize=(6,5));scattermap(ax,d,col,cblue if col.startswith('K32') or col=='I_state_mean' else grey if col=='lower80_mean' else 'viridis',1 if col=='I_state_mean' else limits[col]);ax.set_title(col+'\n'+statuses[sec]);save(fig,o/col,d[['spot_id','x','y',col]],cap)
 fig,axs=plt.subplots(1,4,figsize=(15,3.6))
 for ax,col in zip(axs,['MP10_Top20_UCell','EPCAM_log1pCP10k','nCount_full','Epi_total_mean']):
  ok=d.I_state_mean.notna()&d[col].notna();q=d[ok];ax.scatter(q[col],q.I_state_mean,s=4,alpha=.4,c=blue);rho=float(spearmanr(q[col],q.I_state_mean).statistic) if len(q)>2 and q[col].nunique()>1 else np.nan;ax.set_xlabel(col);ax.set_ylabel('I_state posterior mean');ax.set_ylim(0,1);ax.set_title('Spearman='+format(rho,'.2f'));spcorr.append({'section':sec,'predictor':col,'n_spots':len(q),'spearman':rho,'numerical_status':statuses[sec]})
 save(fig,o/'state_expression_associations',d[['spot_id','I_state_mean','MP10_Top20_UCell','EPCAM_log1pCP10k','nCount_full','Epi_total_mean']],'Per-section descriptive association, no spot-level significance test. Mixed-spot MP10 score is not automatically epithelial intrinsic activity; shared genes are not independent validation.')
wt(pd.DataFrame(clip),F/'display_limits.tsv');wt(pd.DataFrame(spcorr,columns=['section','predictor','n_spots','spearman','numerical_status']),R/'05_evaluation/spatial_associations.tsv')
# Fixed scope support figures.
support=pd.read_csv(R/'00_manifest/support_by_sample_source_group_dataset_split_label_v7.tsv',sep='\t');e=support[support.label_v7.isin([K,L,H])];pv=e.pivot_table(index=['split','dataset','sample'],columns='label_v7',values='n_cells',fill_value=0);fig,ax=plt.subplots(figsize=(12,max(6,len(pv)*.15)));pv[[K,L,H]].plot.barh(stacked=True,ax=ax,color=[blue,'#aaaaaa','#bc6c25']);ax.set_xlabel('Frozen epithelial cells');ax.set_ylabel('split / dataset / sample');ax.tick_params(axis='y',labelsize=5);ax.legend(fontsize=7);save(fig,F/'reference_source_support',e,'Frozen K32 retains low-UCell members. Withheld high non-K32 is omitted only from reference; DEV withheld cells are retained in generalization bags. Conservative source groups are not independently verified patients.')
# Reference diagnostics when available.
ref=R/'01_reference'
if (ref/'signatures.tsv.gz').exists():
 sig=pd.read_csv(ref/'signatures.tsv.gz',sep='\t',index_col=0);old=pd.read_csv(S/'03_reference/MAIN/signatures.tsv.gz',sep='\t',index_col=0);assert sig.index.equals(old.index);corr=np.log1p(sig).corr();wt(corr.reset_index(),ref/'log1p_Pearson_matrix.tsv');wt(np.log1p(sig).corr(method='spearman').reset_index(),ref/'log1p_Spearman_matrix.tsv');fig,ax=plt.subplots(figsize=(12,10));im=ax.imshow(corr,vmin=0,vmax=1,cmap='viridis');ax.set_xticks(range(31),corr.columns,rotation=90,fontsize=6);ax.set_yticks(range(31),corr.index,fontsize=6);fig.colorbar(im,ax=ax);save(fig,F/'reference_log1p_correlation',corr.reset_index(),'Pearson correlations of log1p posterior-mean NB signatures across all11230 G0 genes. Descriptive diagnostic, no fixed correlation pass threshold.')
 effects=pd.read_csv(ref/'TRAIN_expression_effects.tsv.gz',sep='\t').set_index('gene');effects['NB_K32']=sig[K];effects['NB_lower80']=sig[L];effects['NB_log2FC_pc0.1']=np.log2((sig[K]+.1)/(sig[L]+.1));effects['old_K32']=old.Epi_K32;effects['old_other']=old.Epi_other;wt(effects.reset_index(),ref/'gene_signature_comparison.tsv.gz');rows=[]
 for lab in sig.columns:
  ol='Epi_K32' if lab==K else 'Epi_other' if lab==L else lab;rows.append({'new_label':lab,'old_label':ol,'log1p_Pearson':np.corrcoef(np.log1p(sig[lab]),np.log1p(old[ol]))[0,1],'new_sum':sig[lab].sum(),'old_sum':old[ol].sum(),'amplitude_ratio':sig[lab].sum()/old[ol].sum()})
 wt(pd.DataFrame(rows),ref/'old_new_signature_amplitudes.tsv');wt(corr[K].drop(K).sort_values(ascending=False).rename('Pearson_log1p').rename_axis('label').reset_index(),ref/'K32_nearest_signatures.tsv');js(ref/'REFERENCE_DIAGNOSTIC_SUMMARY.json',{'old_K32_other_Pearson_log1p':float(np.corrcoef(np.log1p(old.Epi_K32),np.log1p(old.Epi_other))[0,1]),'new_K32_lower80_Pearson_log1p':float(corr.loc[K,L]),'interpretation':'descriptive; not biological identifiability pass; lower80 construction enriches MP10 contrast'});fig,ax=plt.subplots(figsize=(6,5));ax.scatter(np.log1p(sig[L]),np.log1p(sig[K]),s=2,c='#999999',alpha=.3);disp=pd.read_csv(R/'00_manifest/DISPLAY_GENES.tsv',sep='\t').gene
 for g in disp:ax.scatter(np.log1p(sig.loc[g,L]),np.log1p(sig.loc[g,K]),color=blue);ax.annotate(g,(np.log1p(sig.loc[g,L]),np.log1p(sig.loc[g,K])),fontsize=8)
 ax.set_xlabel('log1p lower80 signature');ax.set_ylabel('log1p K32 anchor signature');save(fig,F/'reference_K32_lower80_scatter',effects.reset_index(),'All G0 genes; highlighted display genes frozen from TRAIN before new spatial predictions. Lower80 selection makes MP10 enrichment partly expected by construction; no new DEG p-values.')
 top=effects[effects.MP10_Top20].sort_values('NB_log2FC_pc0.1',ascending=False);fig,ax=plt.subplots(figsize=(9,5));xx=np.arange(len(top));ax.bar(xx-.2,top.NB_K32,.4,label=K,color=blue);ax.bar(xx+.2,top.NB_lower80,.4,label=L,color='#aaaaaa');ax.set_xticks(xx,top.index,rotation=60,ha='right');ax.set_ylabel('Posterior mean signature');ax.legend(fontsize=7);save(fig,F/'reference_MP10_gene_expression',top.reset_index(),'Frozen MP10 Top20 genes measured in G0; absent G0 genes are not represented as zero. These differences are reference diagnostics, not independent evidence.')
 dr=pd.read_csv(ref/'TRAIN_source_directions.tsv',sep='\t');ss=dr[dr.eligible].groupby('gene').agg(n_supported_samples=('sample','nunique'),n_source_groups=('source_group','nunique'),n_datasets=('dataset','nunique'),positive_samples=('log2FC_pc0.1',lambda x:(x>0).sum()),median_log2FC=('log2FC_pc0.1','median')).reset_index();wt(ss,ref/'TRAIN_source_direction_summary.tsv');pp=dr.pivot(index='sample',columns='gene',values='log2FC_pc0.1');fig,ax=plt.subplots(figsize=(11,7));im=ax.imshow(np.ma.masked_invalid(pp),cmap='coolwarm',vmin=-2,vmax=2,aspect='auto');ax.set_xticks(range(len(pp.columns)),pp.columns,rotation=90,fontsize=7);ax.set_yticks(range(len(pp)),pp.index,fontsize=5);fig.colorbar(im,ax=ax,label='Raw mean log2FC (+0.1)');save(fig,F/'TRAIN_source_gene_directions',dr,'Within TRAIN sample effect directions require >=30 cells in each state. Missing support or absent G0 gene shown NA; descriptive counts, not cell-level significance tests.')
# DEV figures, including each source and dataset.
for model in ['DEV','DEV_OLD']:
 p=R/'05_evaluation'/(model+'_panel_results.tsv.gz')
 if not p.exists():continue
 d=pd.read_csv(p,sep='\t');status=d.numerical_status.iloc[0]
 for ds in ['ALL']+sorted(d.dataset.unique()):
  z=d if ds=='ALL' else d[d.dataset.eq(ds)];o=F/'DEV'/model/ds;ac=z[z.panel.isin(['A','C'])].drop_duplicates('pseudospot_id');fig,axs=plt.subplots(1,3,figsize=(12,3.8))
  for ax,t,pred in zip(axs,['n_K32_true','n_epi_true','K32_fraction_epi_true'],['K32_mean','Epi_total_mean','I_state_mean']):
   for sample,q in ac.groupby('sample'):ax.scatter(q[t],q[pred],s=10,alpha=.65,label=sample)
   hi=max(ac[t].max(),ac[pred].max());ax.plot([0,hi],[0,hi],'--',c='black',lw=.8);ax.set_xlabel('True '+t);ax.set_ylabel(pred)
  axs[-1].legend(fontsize=4,bbox_to_anchor=(1,1));fig.suptitle(ds+' | A/C | '+status);save(fig,o/'quantity_calibration',ac,'A/C unique bags, exact frozen K32 and total epithelial truth. Ratio is posterior draw-level K32/epithelial. Source points are technical repetitions, not patients. '+status)
  matchkeys=['sample','source_group','dataset','repeat','n_epi_design','n_dose_design'];qa=z[z.panel.eq('A')][matchkeys];qb=z[z.panel.eq('B')][matchkeys];commonkeys=qa.merge(qb,on=matchkeys,validate='one_to_one');pairedz=z[z.panel.isin(['A','B'])].merge(commonkeys,on=matchkeys,validate='many_to_one');fig,axs=plt.subplots(2,2,figsize=(11,7))
  for row,ne in enumerate([4,12]):
   for col,y in enumerate(['K32_mean','I_state_mean']):
    ax=axs[row,col]
    for (panel,sample,rep),q in pairedz[pairedz.n_epi_design.eq(ne)].groupby(['panel','sample','repeat']):
     q=q.sort_values('n_dose_design');ax.plot(q.n_dose_design/ne,q[y],color=blue if panel=='A' else '#bd6b26',alpha=.35,marker='o',ms=2,lw=.7)
    ax.set_xlabel('Designed high-pool fraction');ax.set_ylabel(y);ax.set_title('Epithelial='+str(ne));ax.plot([],[],color=blue,label='A: frozen K32');ax.plot([],[],color='#bd6b26',label='B: withheld high');ax.legend(fontsize=7)
  fig.suptitle(ds+' | fixed same-sample backgrounds | '+status);save(fig,o/'A_B_dose_response',pairedz,'Dose series use the same concrete background within sample/repeat/epithelial count. B measures withheld-high state generalization and separately has exact K32 truth0. Only jointly feasible A/B source, repeat, epithelial count and dose conditions are shown; exact background hashes match. '+status)
  fig,axs=plt.subplots(2,4,figsize=(15,7))
  for row,ne in enumerate([4,12]):
   for col,(xx,yy) in enumerate([('high_fraction_epi_true','K32_mean'),('high_fraction_epi_true','I_state_mean'),('mean_ucell_epi_true','K32_mean'),('mean_ucell_epi_true','I_state_mean')]):
    ax=axs[row,col]
    for (panel,sample,rep),qq in pairedz[pairedz.n_epi_design.eq(ne)].groupby(['panel','sample','repeat']):
     qq=qq.sort_values('n_dose_design');ax.plot(qq[xx],qq[yy],color=blue if panel=='A' else '#bd6b26',alpha=.35,marker='o',ms=2,lw=.7)
    ax.set_xlabel(xx);ax.set_ylabel(yy);ax.set_title('Epithelial='+str(ne));ax.plot([],[],color=blue,label='A: K32');ax.plot([],[],color='#bd6b26',label='B: withheld high');ax.legend(fontsize=6)
  fig.suptitle(ds+' | observed frozen-cell state truths | '+status);fig.tight_layout();save(fig,o/'A_B_true_state_response',pairedz,'Predictions against actual bag high-UCell epithelial fraction and mean original epithelial UCell, not designed high-pool fraction. Low-UCell frozen K32 contributes exact K32 but not high-UCell truth. Jointly feasible source/background/dose comparisons only. '+status)
  c=z[z.panel.eq('C')].drop_duplicates('pseudospot_id');fig,ax=plt.subplots(figsize=(9,4));samples=sorted(c['sample'].unique())
  for i,sample in enumerate(samples):
   q=c[c['sample'].eq(sample)];ax.scatter(np.full(len(q),i)+(q.n_epi_true.to_numpy()/20-.5)*.5,q.K32_mean,c=q.n_epi_true,cmap='viridis',vmin=0,vmax=20,s=20)
  ax.set_xticks(range(len(samples)),samples,rotation=75,ha='right',fontsize=6);ax.set_ylabel('K32 anchor posterior mean');fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),ax=ax,label='True epithelial count');ax.set_title(status);save(fig,o/'C_low_background',c,'Pure lower80 epithelial and pure non-epithelial negative-state bags, exact K32=0. Posterior abundance distribution, no binary threshold. '+status)
# Numerical diagnostics are separate from biological evidence.
for name,o in [('NB_LOWER80',R/'01_reference')]+[(n,outmap(n)) for n in ['DEV','DEV_OLD']+config()['sections']]:
 hp=o/'training_history.tsv'
 if not hp.exists():continue
 h=pd.read_csv(hp,sep='\t');fig,ax=plt.subplots(figsize=(7,4));ax.plot(h.epoch,h.loss,color=blue,lw=.7);ax.set_xlabel('Actual completed epoch');ax.set_ylabel('Training ELBO loss');ax.set_title(name);save(fig,F/'numerical'/name,h,'Actual completed epoch loss history; numerical stability uses prespecified loss and abundance windows, not a visual biological judgment.')
 ap=o/'abundance_convergence.tsv'
 if ap.exists():
  a=pd.read_csv(ap,sep='\t');fig,axs=plt.subplots(1,2,figsize=(10,3.8))
  for ax,metric,cut in zip(axs,['median_relative_change','q95_relative_change'],[.05,.20]):
   for lab,q in a.groupby('label',sort=False):ax.plot(q.epoch,q[metric],marker='o',ms=3,label=lab)
   ax.axhline(cut,color='black',linestyle='--',lw=.8,label='Prespecified limit');ax.set_xlabel('Actual completed epoch');ax.set_ylabel(metric);ax.legend(fontsize=6)
  fig.suptitle(name);save(fig,F/'numerical'/(name+'_abundance_stability'),a,'Formal abundance stability: abs(new-old)/(old+0.1) using 256-draw means, checkpoint separation at least2000 actual epochs. All three labels must meet median<0.05 and q95<0.20, together with loss criterion. Flat ELBO alone is insufficient. Longer final comparison gaps are retained in plot data.')
wt(pd.DataFrame(captions),F/'FIGURE_CAPTIONS.tsv');say('FIGURES UPDATED',len(captions))
