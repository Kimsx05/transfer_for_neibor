from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import subprocess
for st in ['Regular','Bold']:font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial:style='+st],text=True))
plt.rcParams.update({'font.family':'Arial','font.size':12,'pdf.fonttype':42,'lines.linewidth':1.2,'axes.linewidth':1.1,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.right':False,'axes.spines.top':False})
COL=dict(zip(DS,['#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#a6cee3','#b15928']));PAL=pd.read_csv(S/'tables/niche_palette.tsv',sep='\t').set_index('accepted_label').color.to_dict();INDEX=[]
def short(s):return str(s).replace('HALLMARK_','').replace('BUDDEN_','Budden ').replace('B1_','').replace('B3_','').replace('_distance',' distance').replace('_SIGNALING','').replace('_RESPONSE',' response').replace('_',' ')
def dsname(s):return s.replace('__UNVERIFIED',' (sample units)')
def export(f,name,keys=(),desc=''):
 p=O/'figures'/name;f.savefig(str(p)+'.pdf',bbox_inches='tight',facecolor='white');f.savefig(str(p)+'.png',bbox_inches='tight',dpi=300,facecolor='white');rec=dict(figure=name,pdf=str(p)+'.pdf',png=str(p)+'.png',description=desc)
 if keys:
  for a in keys:a.set_visible(False)
  f.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight',facecolor='white')
  for a in f.axes:a.set_visible(a in keys)
  f.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight',facecolor='white');rec.update(panel_pdf=str(p)+'_panel_only.pdf',legend_pdf=str(p)+'_legend_only.pdf')
 plt.close(f);INDEX.append(rec)
def dataset_legend(f):
 ax=f.add_axes([.79,.15,.2,.7]);ax.axis('off');ax.legend(handles=[Patch(color=COL[d],label=dsname(d)) for d in DS if d!='GSE285715'],frameon=False,loc='center left',fontsize=11);return ax

def main():
 reg=read('endpoints.tsv.gz');co=read('representative_coordinates.tsv');meta=read('unit_metadata.tsv').set_index('record_id');u=read('unit_endpoints.tsv.gz');ds=read('dataset_endpoint_summary.tsv.gz');disc=read('discussion_relations.tsv');h=read('heterogeneity_configurations.tsv.gz');st=read('stability_summary.tsv');scope=read('search_scope.tsv');ur=combo_id([d for d in DS if d!='GSE285715'])
 # 1 Search volume by module and number of whole datasets, without a millions-of-points plot.
 sc=scope.copy();sc['module']=sc.family.str.split('|').str[0];pv=sc.groupby(['module','n_datasets']).n_rows.sum().unstack().fillna(0);f,ax=plt.subplots(figsize=(9,5));im=ax.pcolormesh(np.log10(pv+1),cmap='viridis');ax.set_xticks(np.arange(len(pv.columns))+.5,pv.columns);ax.set_yticks(np.arange(len(pv.index))+.5,pv.index);ax.set_xlabel('Whole datasets in a combination (n)');ax.set_ylabel('Endpoint module');cb=f.colorbar(im,ax=ax);cb.set_label('log10(number of searched rows + 1)');cb.solids.set_rasterized(False);export(f,'01_search_scope',[cb.ax],'All searched configurations, including single datasets; no Inhouse requirement')
 # 2 All seven-dataset descriptive reference; no false unit inference for missing identities.
 cr=read('combination_associations.tsv.gz');rr=reg[reg.module.str.startswith(('B1','B3'))&reg.signature.eq('full')&reg.canonical_program.isin(CORE)];z=cr[cr.endpoint_id.isin(rr.endpoint_id)&cr.combination.eq('D127')].merge(rr,on='endpoint_id');z['axis']=z.module+' '+z.statistic;mat=z.pivot(index='canonical_program',columns='axis',values='estimate').reindex(index=[p for p in CORE if p in z.canonical_program.unique()]);f,ax=plt.subplots(figsize=(11,8));im=ax.pcolormesh(np.ma.masked_invalid(mat.to_numpy()),cmap='RdBu_r',vmin=-.5,vmax=.5);ax.set_xticks(np.arange(len(mat.columns))+.5,[short(x) for x in mat.columns],rotation=60,ha='right');ax.set_yticks(np.arange(len(mat))+.5,[short(x) for x in mat.index]);ax.invert_yaxis();cb=f.colorbar(im,ax=ax);cb.set_label('Dataset-equal r or partial r (labelled by column)');cb.solids.set_rasterized(False);f.subplots_adjust(left=.4,bottom=.3);export(f,'02_all_dataset_reference',[cb.ax],'Seven-dataset descriptive reference; GSE285715 has section-only records')
 del cr
 # 3 Unit points and dataset means, with between-unit t intervals, never spot-level scatter.
 forests=disc[disc.selection.eq('all_identified_units')&disc.module.isin(['B3_Epi_distance','C1','C2'])].drop_duplicates('module')
 for rr in forests.itertuples():
  z=u[u.endpoint_id.eq(rr.endpoint_id)&u.value.notna()].copy()
  if rr.n_min_rule:z=z[z.n_min.ge(rr.n_min_rule)&z.n_other_min.ge(rr.n_min_rule)]
  f,ax=plt.subplots(figsize=(9,5));rng=np.random.default_rng(SEED)
  for i,d in enumerate(DS):
   w=z[z.dataset_id.eq(d)];v=transform(w.value.to_numpy(),rr.kind);ax.scatter(w.value,i+rng.uniform(-.14,.14,len(w)),s=24,color=COL[d],facecolors='none' if (w.identity_level=='section_only').any() else COL[d],alpha=.65)
   if len(v):
    m=v.mean();se=v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 and w.identity_level.eq('analysis_unit').all() else np.nan;crit=stats.t.ppf(.975,len(v)-1);lo,hi=inverse(np.array([m-crit*se,m+crit*se]),rr.kind);point=inverse(m,rr.kind);ax.plot([lo,hi],[i,i],color='black',lw=1.3);ax.plot(point,i,'D',color='black',ms=6)
  ax.set_yticks(range(len(DS)),[dsname(d)+(' (section-only)' if d=='GSE285715' else '') for d in DS]);ax.axvline(0,color='grey',ls='--',lw=1);ax.set_xlabel(('Paired delta r' if rr.module=='C2' else 'Unit r')+'\n'+short(rr.canonical_program)+((' / '+rr.candidate_type) if rr.candidate_type!='none' else ''));f.subplots_adjust(left=.33,bottom=.2);export(f,'03_forest_'+rr.module,desc='Points are actual unit effects or explicitly section-only records; black intervals use between-record uncertainty')
 # 4 Compare explained variance, dataset dependence and practical cluster structure.
 names=co.input_family.unique();best=[]
 for name in names:
  z=h[h.input_family.eq(name)&h.combination.eq(ur)].sort_values('silhouette',ascending=False);best.append(z.iloc[0])
 best=pd.DataFrame(best);f,ax=plt.subplots(figsize=(10,5));x=np.arange(len(best));ax.bar(x-.22,best.PC1_variance,.22,label='PC1 variance',color='#1f78b4');ax.bar(x,best.PC2_variance,.22,label='PC2 variance',color='#a6cee3');ax.bar(x+.22,best.PC1_dataset_eta2,.22,label='Dataset eta² on PC1',color='#b15928');ax.set_xticks(x,[short(n) for n in best.input_family],rotation=30,ha='right');ax.set_ylabel('Variance fraction');l=f.add_axes([.25,.02,.7,.06]);l.axis('off');l.legend(handles=[Patch(color=c,label=v) for c,v in zip(['#1f78b4','#a6cee3','#b15928'],['PC1 variance','PC2 variance','Dataset eta² on PC1'])],ncol=3,frameon=False);f.subplots_adjust(bottom=.38);export(f,'04_input_heterogeneity',[l])
 # 5 Five fixed input representatives, PC1/PC2 only; dataset colors are consistent.
 for i,name in enumerate(names):
  z=co[co.input_family.eq(name)].merge(meta.reset_index(),on='record_id');b=best[best.input_family.eq(name)].iloc[0];f,ax=plt.subplots(figsize=(10,5.5));
  for d,w in z.groupby('dataset_id'):ax.scatter(w.PC1,w.PC2,color=COL[d],s=35,edgecolors='white',linewidths=.5)
  ax.set_xlabel(f'PC1 ({100*b.PC1_variance:.1f}% variance)');ax.set_ylabel(f'PC2 ({100*b.PC2_variance:.1f}% variance)');f.subplots_adjust(right=.76);l=dataset_legend(f);export(f,f'05_PCA_{i+1}',[l],name+'; no class claim from PCA separation')
 # 6 Mean labelled-Epi niche proportions in the internally derived Ward candidate groups.
 bg=read('unit_background.tsv.gz');bg=bg[bg.input_family.eq('niche_all_Epi_all')];labels=co[co.input_family.eq('niche_all_Epi_all')][['record_id','cluster']];z=bg.merge(labels,on='record_id');mat=z.groupby(['cluster','feature']).value.mean().unstack().reindex(columns=[str(i) for i in range(9)]);f,ax=plt.subplots(figsize=(8,5));base=np.zeros(len(mat));
 for k in range(9):ax.bar(np.arange(len(mat)),mat[str(k)]*100,bottom=base*100,color=PAL[k],width=.7,label=str(k));base+=mat[str(k)].to_numpy()
 ax.set_xticks(np.arange(len(mat)),['Candidate '+str(int(k)+1)+f' (n={(labels.cluster==k).sum()})' for k in mat.index],rotation=15,ha='right');ax.set_ylabel('Mean labelled-Epi niche composition (%)');l=f.add_axes([.82,.15,.15,.7]);l.axis('off');l.legend(handles=[Patch(color=PAL[k],label=str(k)) for k in range(9)],title='Formal R12 label',frameon=False);f.subplots_adjust(right=.8,bottom=.22);export(f,'06_niche_candidate_composition',[l],'Groups are described by their construction variables, not independently validated')
 # 7 Stability, preserving the rare 53+1 R0 solution rather than presenting it as two robust types.
 ks=['unit_subset_80percent','feature_subset_80percent','kmeans_same_K'];mat=st[st.check.isin(ks)].pivot(index='input_family',columns='check',values='mean').reindex(index=names,columns=ks);f,ax=plt.subplots(figsize=(10,5));im=ax.pcolormesh(mat,cmap='viridis',vmin=0,vmax=1);ax.set_yticks(np.arange(len(mat))+.5,[short(x) for x in mat.index]);ax.set_xticks(np.arange(3)+.5,['80% unit subset','80% feature subset','K-means sensitivity'],rotation=20,ha='right');cb=f.colorbar(im,ax=ax);cb.set_label('Mean adjusted Rand index');cb.solids.set_rasterized(False);f.subplots_adjust(left=.4,bottom=.23);export(f,'07_internal_stability',[cb.ax])
 # 8 Actual unit residual plots for the strongest nominal partial ecological relation and a paired-delta example.
 from ecology import design
 eco=read('ecology_associations.tsv.gz').merge(reg,on='endpoint_id');pc=co[co.input_family.eq('niche_all_Epi_all')].set_index('record_id')[['PC1','PC2']];pc=(pc-pc.mean())/pc.std();specs=[('C1','vCAF','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','partial_r'),('C2','NK','HALLMARK_DNA_REPAIR','paired_delta_r')]
 for j,(mo,ca,pr,stt) in enumerate(specs):
  er=eco[eco.module.eq(mo)&eco.candidate_type.eq(ca)&eco.canonical_program.eq(pr)&eco.statistic.eq(stt)&eco.score_method.eq('AddModuleScore_SCT_data')&eco.signature.eq('full')&eco.ecology_axis.eq('PC1')].iloc[0];v=u[u.endpoint_id.eq(er.endpoint_id)&u.identity_level.eq('analysis_unit')&u.n_min.ge(10)&u.n_other_min.ge(10)&u.value.notna()&u.record_id.isin(pc.index)];xx=pc.loc[v.record_id];yy=transform(v.value.to_numpy(),er.kind);dd=v.dataset_id.to_numpy();cov=design(xx[['PC2']].to_numpy(),dd);xr=xx.PC1.to_numpy()-cov@np.linalg.lstsq(cov,xx.PC1.to_numpy(),rcond=None)[0];yr=yy-cov@np.linalg.lstsq(cov,yy,rcond=None)[0];f,ax=plt.subplots(figsize=(10,5.5));
  for d in sorted(set(dd)):ok=dd==d;ax.scatter(xr[ok],yr[ok],c=COL[d],s=35)
  xs=np.array([xr.min(),xr.max()]);ax.plot(xs,np.polyval(np.polyfit(xr,yr,1),xs),color='black');ax.set_xlabel('Niche ecological PC1 (residual within dataset, controlling PC2)');ax.set_ylabel('Residual '+('Fisher z (partial r)' if mo=='C1' else 'paired delta r')+'\n'+ca+' / '+short(pr));f.subplots_adjust(left=.18,right=.76,bottom=.18);l=dataset_legend(f);export(f,f'08_ecology_{j+1}',[l],'Unit-level residuals; nominal association did not pass full ecological-search BH')
 # 9 Raw/partial and full/noMP10 sensitivity of B3 programme gradients.
 se=read('discussion_sensitivity.tsv.gz');se=se[se.module.str.startswith('B3')&se.combination.eq(ur)];lab=se.module+' / '+se.actual_program+' / '+se.statistic;f,ax=plt.subplots(figsize=(10,max(5,.32*len(se))));y=np.arange(len(se));ax.errorbar(se.estimate,y,xerr=np.vstack([se.estimate-se.conditional_low,se.conditional_high-se.estimate]),fmt='o',color='#1f78b4',ecolor='#777777',capsize=2);ax.set_yticks(y,[short(x) for x in lab]);ax.axvline(0,color='grey',ls='--');ax.set_xlabel('Dataset-equal r / partial r (unit-level conditional interval)');ax.invert_yaxis();f.subplots_adjust(left=.62);export(f,'09_B3_sensitivity',desc='All predesignated variants for the selected programmes, not a method picked by effect')
 save(INDEX,'figure_index.tsv');print('figures complete',len(INDEX),flush=True)
if __name__=='__main__':main()
