from fu_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
fm.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1,'lines.linewidth':1.2,'pdf.fonttype':42,'figure.facecolor':'white','savefig.facecolor':'white'})
spec=importlib.util.spec_from_file_location('vp',V/'code/10_plot.py');vp=importlib.util.module_from_spec(spec);spec.loader.exec_module(vp);vp.FIG=R/'figures';save=vp.save
read=lambda f:pd.read_csv(R/'tables'/f,sep='\t')
# Post-hoc truth attribution only; no expression, score, threshold, or evidence change.
m=pd.read_parquet(R/'tables/expanded_mixtures.parquet');z=pd.read_parquet(R/'tables/mixture_rule_results.parquet');assert np.array_equal(m.Treg_main,z.Treg_main)
for c in ['realized_T_RNA_UMI','realized_Treg_RNA_UMI','realized_Treg_total_RNA_fraction','realized_Treg_T_RNA_fraction','truth_definition']:z[c]=m[c].to_numpy()
z.to_parquet(R/'tables/mixture_rule_results.parquet',index=False)
d=read('paired_effects.tsv');fig,axs=plt.subplots(1,2,figsize=(12,5));names=[]
for i,r in enumerate(d.itertuples()):
 names.append(r.variant.replace('_',' ')+' / '+str(r.block_um)+' µm')
 for col,off,key in [('#397da5',-.09,'beta_bin'),('#c27a44',.09,'beta_cell')]:axs[0].plot([getattr(r,key+'_low'),getattr(r,key+'_high')],[i+off]*2,c=col);axs[0].plot(getattr(r,key),i+off,'o',c=col)
 axs[1].plot([r.delta_low,r.delta_high],[i]*2,c='#6e568a');axs[1].plot(r.delta,i,'o',c='#6e568a')
for ax in axs:ax.axvline(0,c='gray',ls='--');ax.set_yticks(range(len(names)),names);ax.invert_yaxis()
axs[0].set_xlabel('Standardized rank β (paired block 95% CI)');axs[1].set_xlabel('Δβ = cell − bin (paired block 95% CI)');axs[1].set_yticklabels([]);hh=[Line2D([],[],c='#397da5',marker='o',label='Bin source'),Line2D([],[],c='#c27a44',marker='o',label='Cell source')];axs[0].legend(handles=hh,frameon=False,loc='lower left',bbox_to_anchor=(0,1.01),ncol=2);fig.tight_layout();save(fig,'Cervilla_paired_source_effects',legends=hh)
d=read('stepwise_decomposition.tsv');fig,ax=plt.subplots(figsize=(7,4));names=[]
for i,r in enumerate(d.itertuples()):
 names.append(r.step.replace('_',' ')+' / '+r.unit);ax.plot([r.low,r.high],[i]*2,c='#397da5');ax.plot(r.effect,i,'o',c='#397da5')
ax.set_yticks(range(len(d)),names);ax.invert_yaxis();ax.axvline(0,c='gray',ls='--');ax.set_xlabel('Standardized rank β (95% block CI)');save(fig,'Cervilla_stepwise_decomposition')
d=read('leave_one_block_out.tsv');fig,axs=plt.subplots(1,2,figsize=(10,4));hh=[]
for ax,rad in zip(axs,[400,800]):
 q=d[d.block_um.eq(rad)]
 for col,key in [('#397da5','beta_bin'),('#c27a44','beta_cell'),('#6e568a','delta')]:
  h=ax.plot(range(len(q)),q[key],'o-',ms=3,c=col,label=key)[0]
  if rad==400:hh.append(h)
 ax.axhline(0,c='gray',ls='--');ax.set_xlabel(f'Omitted spatial block ({rad} µm)');ax.set_ylabel('Fixed-design rank effect')
axs[1].legend(handles=hh,frameon=False);save(fig,'Cervilla_leave_one_block_out',legends=hh)
d=read('source_coverage_summary.tsv');a=pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet');cov=read('Cervilla_source_coverage.tsv');cov=cov[cov.radius_um.eq(200)&cov.T_region_id.isin(a.T_region_id)];fig,axs=plt.subplots(1,3,figsize=(12,3.8));axs[0].scatter(a.MP10_bin,a.MP10_cell,s=10,c='#397da5');axs[0].set_xlabel('Bin MP10 mean');axs[0].set_ylabel('Cell MP10 mean');axs[1].hist(cov.Jaccard,bins=np.linspace(0,1,21),color='#6e568a');axs[1].set_xlabel('8 µm source-footprint Jaccard');axs[1].set_ylabel('Common T centers');q=read('common_support_paired_effects.tsv')
for i,r in enumerate(q.itertuples()):axs[2].plot([r.delta_low,r.delta_high],[i]*2,c='#6e568a');axs[2].plot(r.delta,i,'o',c='#6e568a')
axs[2].set_yticks([0,1],['400 µm blocks','800 µm blocks']);axs[2].axvline(0,c='gray',ls='--');axs[2].set_xlabel('Reconstructed support Δβ (n=170)');fig.tight_layout();save(fig,'Cervilla_common_support')
# All five unchanged primary results, including descriptive-only branches.
d=read('v1_five_section_main.tsv');fig,ax=plt.subplots(figsize=(8,5));hh=[]
for radius,col,off in [(50,'#397da5',-.12),(200,'#c27a44',.12)]:
 hh.append(Line2D([],[],c=col,marker='o',label=f'{radius} µm (v1 unchanged)'))
 for i,s in enumerate(SECTIONS):
  r=d[d.section.eq(s)&d.radius_um.eq(radius)].iloc[0];ax.plot(r.effect,i+off,'o',c=col,mfc=col if np.isfinite(r.low)else'white')
  if np.isfinite(r.low):ax.plot([r.low,r.high],[i+off]*2,c=col)
ax.axvline(0,c='gray',ls='--');ax.set_yticks(range(5),[s+(' (mask unknown)'if s.startswith('MPBC')else'')for s in SECTIONS]);ax.invert_yaxis();ax.set_xlabel('v1 adjusted rank effect (95% block CI)');hh.append(Line2D([],[],c='gray',marker='o',mfc='white',ls='',label='Descriptive (no reliable CI/P)'));ax.legend(handles=hh,loc='upper left',bbox_to_anchor=(1,1),frameon=False);save(fig,'five_section_v1_forest',legends=hh)
flow=read('qualified_HD_flow.tsv');fig,ax=plt.subplots(figsize=(8,4.5));hh=[]
for i,s in enumerate(SECTIONS):
 q=flow[flow.section.eq(s)];v=np.r_[q.qualified_n.iloc[0],q.survives];h=ax.plot(range(5),100*v/v[0],'o-',label=s)[0];hh.append(h)
ax.set_xticks(range(5),['T eligible','≥3 main','≥1 core','Score >0','Old high'],rotation=20);ax.set_ylabel('Percent of eligible T32 regions');ax.set_yscale('log');ax.legend(handles=hh,loc='upper left',bbox_to_anchor=(1,1),frameon=False);save(fig,'five_section_joint_rule_flow',legends=hh)
perf=read('mixture_per_condition_performance.tsv');p=perf[perf.partition.eq('sample_holdout')&perf.metric.isin(['joint_old','joint_candidate'])];fig,axs=plt.subplots(1,3,figsize=(14,4),sharey=True);hh=[]
for ax,dep in zip(axs,[500,2000,5000]):
 for j,bg in enumerate(['Epi40','Epi70','Epi90','Epi40_Fib30','Epi40_Myeloid30']):
  for metric,ls in [('joint_old','-'),('joint_candidate','--')]:
   q=p[p.depth.eq(dep)&p.background.eq(bg)&p.metric.eq(metric)].groupby('Treg_T_RNA_fraction').rate.mean();ax.plot(q.index,q.values,c=plt.cm.tab10(j),ls=ls,marker='o',ms=3)
  if dep==500:hh.append(Line2D([],[],c=plt.cm.tab10(j),label=bg))
 ax.set_xlabel(f'Treg / T RNA fraction ({dep} total UMI)');ax.set_ylim(-.01,1)
axs[0].set_ylabel('Held-out full-rule detection (3 competitors mean)');hh.extend([Line2D([],[],c='black',ls='-',label='Old'),Line2D([],[],c='black',ls='--',label='Candidate (not adopted)')]);axs[-1].legend(handles=hh,loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=10);save(fig,'mixture_detectability',legends=hh)
# Background-specific negative controls; no pooling conceals Tfh or activation.
p=perf[perf.partition.eq('sample_holdout')&perf.Treg_T_RNA_fraction.eq(0)&perf.metric.isin(['joint_old','joint_candidate'])];fig,axs=plt.subplots(1,3,figsize=(14,4),sharey=True);hh=[]
for ax,comp in zip(axs,p.competitor.unique()):
 for j,metric in enumerate(['joint_old','joint_candidate']):
  q=p[p.competitor.eq(comp)&p.metric.eq(metric)].sort_values(['background','depth']);x=np.arange(len(q))+(j-.5)*.16;ax.plot(x,q.rate,'o',c=['#397da5','#c27a44'][j],label=metric);ax.vlines(x,q.low,q.high,colors=['#397da5','#c27a44'][j],lw=1)
 ax.axhline(.05,c='gray',ls='--');ax.set_xticks(range(len(q)),[f'{r.background}\n{r.depth}'for r in q.itertuples()],rotation=90,fontsize=7);ax.set_xlabel(comp,fontsize=10)
axs[0].set_ylabel('Held-out full-rule FPR (conditional Wilson CI)');hh=[Line2D([],[],c=c,marker='o',ls='',label=x)for c,x in [('#397da5','Old'),('#c27a44','Candidate (not adopted)')]];axs[-1].legend(handles=hh,frameon=False,fontsize=9);save(fig,'mixture_negative_controls',legends=hh)
# Five whole-section evidence diagnostics. Blank is no source coverage, not zero.
for s in SECTIONS:
 t=pd.read_parquet(R/'tables'/f'{s}_evidence.parquet');fig,axs=plt.subplots(1,2,figsize=(13,5));good=t.measurement_eligible;zero=good&t.Treg_main.eq(0)&t.Treg_main_detected.gt(0);pos=good&t.Treg_main.gt(0);old=t.old_evidence.eq('Treg_like_region')
 # Use exact existing evidence string, recorded below, rather than changing identity.
 old=t.old_evidence.str.contains('Treg',case=False)&~t.old_evidence.str.contains('insufficient',case=False)
 colors=[('#dddddd',~good,'Unmeasurable T32'),('#7893ab',good,'T measurement eligible'),('#be5a40',old,'Old Treg-like (unchanged)')]
 hh=[]
 for col,k,label in colors:axs[0].scatter(t.loc[k,'x_um'],t.loc[k,'y_um'],s=1.2 if label!='Old Treg-like (unchanged)' else 9,c=col,lw=0,rasterized=True);hh.append(Patch(color=col,label=label))
 for col,k,label in [('#dddddd',~good,'Unmeasurable T32'),('#735b9b',zero,'Detected signature, score zero'),('#45919c',pos,'Eligible positive score')]:axs[1].scatter(t.loc[k,'x_um'],t.loc[k,'y_um'],s=1.3,c=col,lw=0,rasterized=True);hh.append(Patch(color=col,label=label))
 for ax in axs:ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('Relative x (µm); mask unknown'if s.startswith('MPBC')else'x (µm)');ax.set_ylabel('y (µm)')
 axs[-1].legend(handles=hh,loc='upper left',bbox_to_anchor=(1,1),frameon=False,fontsize=9);save(fig,s+'_evidence_diagnostics',legends=hh)
print('SUMMARY_PLOTS_DONE',flush=True)
