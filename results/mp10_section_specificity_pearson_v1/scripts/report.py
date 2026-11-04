from aggregate import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from PIL import Image
import subprocess,gzip,json,textwrap
for style in ['Regular','Bold']:font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial:style='+style],text=True))
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.titlesize':12,'axes.labelsize':12,'pdf.fonttype':42,'lines.linewidth':1.2,'axes.linewidth':1.1,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.right':False,'axes.spines.top':False})
PAL=pd.read_csv(O/'tables/niche_palette.tsv',sep='\t').set_index('accepted_label').color.to_dict();COL=['#D73027','#1F78B4'];RINGS=['R0','R1','R2','R3_R5','R6_R10']
def short(s):return str(s).replace('HALLMARK_','').replace('BUDDEN_','Budden ').replace('EPITHELIAL_MESENCHYMAL_TRANSITION','EMT').replace('FLORES_CRISPR_PRO_INVASION_10','Flores pro-invasion 10').replace('PURAM_PEMT_TOP30','Puram p-EMT Top30').replace('__NO_MP10',' (no MP10)').replace('_',' ')
def export(fig,name,keyaxes=()):
 p=O/'figures'/name;p.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(p)+'.pdf',facecolor='white',bbox_inches='tight');fig.savefig(str(p)+'.png',dpi=300,facecolor='white',bbox_inches='tight')
 if keyaxes:
  for a in keyaxes:a.set_visible(False)
  fig.savefig(str(p)+'_panel_only.pdf',facecolor='white',bbox_inches='tight')
  for a in fig.axes:a.set_visible(a in keyaxes)
  for t in fig.texts:t.set_visible(False)
  fig.savefig(str(p)+'_legend_only.pdf',facecolor='white',bbox_inches='tight')
 plt.close(fig);return str(p)+'.pdf'
def hm(ax,df,title,lim=1):
 cm=plt.get_cmap('RdBu_r').copy();cm.set_bad('#dddddd');im=ax.pcolormesh(np.ma.masked_invalid(df.to_numpy(float)),cmap=cm,vmin=-lim,vmax=lim);ax.set_xticks(np.arange(len(df.columns))+.5,[short(c) for c in df.columns],rotation=60,ha='right');ax.set_yticks(np.arange(len(df.index))+.5,[short(p) for p in df.index]);ax.invert_yaxis();ax.set_title(title);return im

def selected(de,a,normal=False):
 if normal:
  v=a[a.population.eq('normal_Epi')&a.range_type.eq('full_range')&a.score_method.eq('AddModuleScore_SCT_data')&a.program.isin(CORE)&a.n_spots.ge(30)].copy();v['priority']=v.pearson_r.abs()
 else:
  v=de[de.range_type.eq('common_support')&de.score_method.eq('AddModuleScore_SCT_data')&de.program.isin(CORE)&de.n_P95.ge(30)&de.n_nonP95.ge(30)&de.n_blocks_P95.ge(8)&de.n_blocks_nonP95.ge(8)].copy();v['priority']=v.delta_r.abs()
 v=v.sort_values('priority',ascending=False);out=[]
 if len(v):out.append((v.iloc[0].candidate_type,v.iloc[0].program))
 if not normal and len(v):
  opposite=v[v.delta_r*v.iloc[0].delta_r<0]
  if len(opposite):out.append((opposite.iloc[0].candidate_type,opposite.iloc[0].program))
 for r in v.itertuples():
  if len(out)>=2:break
  if (r.candidate_type,r.program) not in out:out.append((r.candidate_type,r.program))
 for pair in [('Treg','BUDDEN_S1P_RESPONSE_18'),('Matrix_CAF','BUDDEN_TNF_EMT_RESPONSE_15')]:
  if len(out)<2 and pair not in out:out.append(pair)
 return out

def scatter_pair(axs,d,a,de,t,p,normal=False,range_type="common_support"):
 x=d['d_candidate__'+t].to_numpy();y=d['AMS_SCT__'+p].to_numpy();uc=d[p].to_numpy();valid=np.isfinite(x)&np.isfinite(y)&np.isfinite(uc);h=d.MP10_P95.to_numpy();low=d.tumor_epi.to_numpy()&~h
 rang=range_type;dd=de[de.candidate_type.eq(t)&de.program.eq(p)&de.score_method.eq('AddModuleScore_SCT_data')&de.range_type.eq(rang)] if len(de) else pd.DataFrame()
 if normal:groups=[('normal_Epi',d['query'].to_numpy()),('normal_Epi',d['query'].to_numpy())];rang='full_range'
 else:groups=[('P95',h),('nonP95',low)]
 if rang!='full_range' and len(dd) and np.isfinite(dd.iloc[0].common_min_um) and np.isfinite(dd.iloc[0].common_max_um):valid&=(x>=dd.iloc[0].common_min_um)&(x<=dd.iloc[0].common_max_um)
 else:
  if not normal and rang!='full_range':valid[:]=False
 allok=valid&(groups[0][1]|groups[1][1]);xr=(x[allok].min(),x[allok].max()) if allok.any() else (0,1000);yr=(y[allok].min(),y[allok].max()) if allok.any() else (-.2,.2);pad=max(.01,(yr[1]-yr[0])*.1)
 for j,((pop,mask),ax) in enumerate(zip(groups,axs)):
  ok=valid&mask;ax.scatter(x[ok],y[ok],s=14,c=np.where(x[ok]==0,'#ff9900',COL[j]),alpha=.65,linewidths=0)
  if ok.sum()>=2 and np.std(x[ok])>0:
   fit=np.polyfit(x[ok],y[ok],1);xx=np.array(xr);ax.plot(xx,np.polyval(fit,xx),color='black',lw=1.4)
  rr=a[a.candidate_type.eq(t)&a.program.eq(p)&a.population.eq(pop)&a.range_type.eq(rang)&a.score_method.eq('AddModuleScore_SCT_data')]
  info='Not estimable'
  if len(rr):
   r=rr.iloc[0];info=f'r={r.pearson_r:.3f}; P={r.p_pearson:.2g}; n={int(r.n_spots)}\nSpatial 95% CI [{r.spatial_ci_low:.3f}, {r.spatial_ci_high:.3f}]'
  ax.set_title(f'{pop} | {t}\n{short(p)}');ax.text(.02,.98,info,transform=ax.transAxes,va='top',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'));ax.set_xlim(xr[0]-5,xr[1]+5);ax.set_ylim(yr[0]-pad,yr[1]+pad);ax.set_xlabel('Nearest candidate P95 core (μm)');ax.set_ylabel('AddModuleScore (SCT data)')
  if normal and j==1:
   ax.clear();ax.text(.1,.6,'Adjacent normal reference\nTumor MP10-P95 is not defined.',transform=ax.transAxes);ax.axis('off')
 if len(dd):
  r=dd.iloc[0]
  if rang=='full_range':return f'{t} / {short(p)}: descriptive full-range delta r={r.delta_r:.3f}; each group retains its full observed range.'
  return f'{t} / {short(p)}: delta r={r.delta_r:.3f} [{r.spatial_ci_low:.3f}, {r.spatial_ci_high:.3f}]; common {r.common_min_um:.0f}–{r.common_max_um:.0f} μm; spatial P={r.p_spatial_fisher_z:.2g}'
 return 'No tumor common-support contrast; no range is fabricated.'

def report(sec):
 if (O/'reports'/f'{sec}_figure_index.tsv').exists():return sec+' cached'
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');a=pd.read_parquet(O/'tables'/f'{sec}_pearson.parquet');de=pd.read_parquet(O/'tables'/f'{sec}_delta.parquet');r=M.set_index('section_id').loc[sec];normal=r.source_context=='adjacent_normal_source';base='sections/'+r.dataset_id+'/'+sec;pages=[];q=d['query'];tum=d.tumor_epi;hi=d.MP10_P95
 old=pd.read_csv(OLD/'tables'/f'{sec}_effects.tsv.gz',sep='\t');new=pd.read_csv(O/'tables'/f'{sec}_new5_context.tsv.gz',sep='\t');context=pd.concat([old,new],ignore_index=True)
 # Page 1: spatial identity, original context and front-only niche denominators.
 f,axes=plt.subplots(2,3,figsize=(22,14));ax=axes.ravel();keys=[];f.suptitle(sec+' | '+str(r.source_context)+'\nFrozen mixed-spot MP10 / front / R12',fontsize=15)
 for i,(col,title,mask,lim,cmap) in enumerate([(RAW,'MP10 AUCell and frozen P95',q,(0,json.loads((OLD/'tables/plot_scales.json').read_text())[RAW]),'viridis'),('d_front','Strict v3.4 signed distance (μm)',d.d_front.notna(),(-1000,1000),'RdBu_r')]):
  ax[i].scatter(d.x,d.y,s=5,c='#ddd',linewidths=0);ok=np.asarray(mask)&d[col].notna();im=ax[i].scatter(d.loc[ok,'x'],d.loc[ok,'y'],c=d.loc[ok,col],s=6,cmap=cmap,vmin=lim[0],vmax=lim[1],linewidths=0);cb=f.colorbar(im,ax=ax[i],shrink=.7);cb.solids.set_rasterized(False);keys.append(cb.ax);ax[i].set_title(title);ax[i].set_aspect('equal');ax[i].invert_yaxis();ax[i].axis('off')
  if i==0:ax[i].scatter(d.loc[hi,'x'],d.loc[hi,'y'],s=18,facecolors='none',edgecolors='#d73027',linewidths=1)
 ax[2].scatter(d.x,d.y,c=[PAL.get(v,'#ddd') for v in d.accepted_label],s=6,linewidths=0);ax[2].set_aspect('equal');ax[2].invert_yaxis();ax[2].axis('off');ax[2].set_title('Formal R12 0–8; grey = unlabelled')
 rings=context[context.program.isin(CORE)&context.analysis.isin(['B1_'+z for z in RINGS])&context.adjustment.eq('adjusted')];mat=rings.pivot_table(index='program',columns='analysis',values='effect').reindex(index=CORE,columns=['B1_'+z for z in RINGS]);mat.columns=RINGS;im=hm(ax[3],mat,'UCell neighborhood: adjusted MP10 coefficient',.6);cb=f.colorbar(im,ax=ax[3],shrink=.8);keys.append(cb.ax);cb.solids.set_rasterized(False)
 for pop,mask,color in [('All internal Epi',tum,COL[1]),('MP10-P95',hi,COL[0])]:
  xx=d.loc[mask&d.d_front.notna(),['d_front',RAW]].copy();xx['bin']=np.floor(xx.d_front/100)*100;vv=xx.groupby('bin').agg(x=('d_front','median'),y=(RAW,'median'),n=(RAW,'count'));vv=vv[vv.n>=5];ax[4].plot(vv.x,vv.y,'o-',color=color,label=pop)
 ax[4].axvspan(-200,0,color='#ffd180',alpha=.25);ax[4].set_xlabel('Strict signed front distance (μm)');ax[4].set_ylabel('MP10 AUCell median (bins n≥5)');ax[4].set_title('Observed coverage only; no extrapolation')
 # show actual gaps by plotted points; line segments only short adjacent bins
 for line in ax[4].lines:
  xx=line.get_xdata();yy=line.get_ydata();line.set_linestyle('None')
  for j in range(len(xx)-1):
   if xx[j+1]-xx[j]<=180:ax[4].plot(xx[j:j+2],yy[j:j+2],color=line.get_color(),lw=1.2)
 niche=pd.read_csv(O/'tables'/f'{sec}_niche_counts.tsv.gz',sep='\t');z=niche[niche.scope.eq('front_only')]
 for j,pop in enumerate(['P95','nonP95']):
  v=z[z.population.eq(pop)].set_index('niche').reindex(range(9));ax[5].bar(np.arange(9)+(j-.5)*.35,v.fraction*100,width=.35,color=COL[j],label=pop)
 ax[5].set_xticks(range(9));ax[5].set_xlabel('R12 formal niche');ax[5].set_ylabel('Front-only labelled Epi (%)');ax[5].set_title('Frozen P95 vs nonP95 within [-200, 0] μm')
 leg=f.add_axes([.36,.012,.61,.035]);leg.axis('off');leg.legend(handles=[Patch(color=PAL[k],label=str(k)) for k in range(9)]+[Patch(color=COL[0],label='P95'),Patch(color=COL[1],label='nonP95 / all Epi')],ncol=11,frameon=False);keys.append(leg)
 f.text(.015,.018,f'Tissue {len(d)} | Epi {q.sum()} | internal {tum.sum()} | P95 {hi.sum()}\nExact front in Epi {(q&d.d_front.notna()).sum()} | R12 labelled {(q&d.accepted_label.notna()).sum()}',fontsize=12);f.subplots_adjust(left=.18,right=.97,top=.9,bottom=.16,hspace=.42,wspace=.6);pages.append(export(f,base+'_1_context',keys))
 # Page 2: fixed panel for every section, never significance filtered.
 f,axs=plt.subplots(2,1,figsize=(23,17));keys=[];pop='normal_Epi' if normal else 'P95';v=a[a.score_method.eq('AddModuleScore_SCT_data')&a.population.eq(pop)&a.range_type.eq('full_range')];mat=v.pivot(index='program',columns='candidate_type',values='pearson_r').reindex(index=CORE,columns=CAND);im=hm(axs[0],mat,pop+' raw Pearson r | full observed range',1);cb=f.colorbar(im,ax=axs[0],shrink=.8);cb.set_label('Pearson r');cb.solids.set_rasterized(False);keys.append(cb.ax)
 if len(de):v=de[de.score_method.eq('AddModuleScore_SCT_data')&de.range_type.eq('common_support')];mat=v.pivot(index='program',columns='candidate_type',values='delta_r').reindex(index=CORE,columns=CAND)
 else:mat=pd.DataFrame(np.nan,index=CORE,columns=CAND)
 im=hm(axs[1],mat,'P95 minus nonP95 r | common support | grey = not estimable',1);cb=f.colorbar(im,ax=axs[1],shrink=.8);cb.set_label('Delta r (display clipped ±1; values unchanged)');cb.solids.set_rasterized(False);keys.append(cb.ax);f.suptitle(sec+' | all 29 candidates; fixed core programs',fontsize=15);f.subplots_adjust(left=.24,right=.96,top=.94,bottom=.16,hspace=.75);pages.append(export(f,base+'_2_pearson',keys))
 # Page 3: reproducible effect/coverage selection, opposite-direction contrast when available.
 pairs=selected(de,a,normal);f,axs=plt.subplots(2,2,figsize=(18,13));notes=[]
 for i,(t,p) in enumerate(pairs):notes.append(scatter_pair(axs[i],d,a,de,t,p,normal))
 f.suptitle(sec+' | exploratory scatter display selected by coverage and |delta r|, not P',fontsize=14);f.text(.06,.025,'\n'.join(notes)+'\nOrange points: distance = 0. Pearson P assumes independent spots; spatial intervals use 500 μm blocks.',fontsize=11);leg=f.add_axes([.4,.002,.3,.02]);leg.axis('off');leg.legend(handles=[Patch(color=COL[0],label='P95'),Patch(color=COL[1],label='nonP95'),Patch(color='#ff9900',label='distance 0')],ncol=3,frameon=False);f.subplots_adjust(top=.91,bottom=.17,hspace=.36,wspace=.25);pages.append(export(f,base+'_3_scatter',[leg]))
 # R qpdf keeps all original vector elements.
 out=O/'reports'/f'{sec}.pdf';expr='qpdf::pdf_combine(c('+','.join(json.dumps(s) for s in pages)+'),output='+json.dumps(str(out))+')';subprocess.run(['Rscript','-e',expr],check=True,stdout=subprocess.DEVNULL)
 ix=[dict(section_id=sec,dataset_id=r.dataset_id,page=i+1,pdf=p,png=p[:-4]+'.png') for i,p in enumerate(pages)];save(ix,O/'reports'/f'{sec}_figure_index.tsv');save([dict(section_id=sec,candidate_type=t,program=p,selection='fixed_core_coverage_then_absolute_delta_with_opposite_direction_if_available') for t,p in pairs],O/'reports'/f'{sec}_scatter_selection.tsv')
 text=f'# {sec}\n\n来源: {r.source_context}; analysis unit: {r.analysis_unit_uid}; patient: {r.patient_uid}。\n\n有效 Epi {q.sum()}，内部 Epi {tum.sum()}，冻结 P95 {hi.sum()}；有精确前沿距离的 Epi {(q&d.d_front.notna()).sum()}；R12 标注 {(q&d.accepted_label.notna()).sum()}。\n\n全部29候选×85程序均有结果记录，原始 Pearson、普通 P、空间区间分别保存。共同范围之外不作组间正式比较。\n\n'+ '\n\n'.join(notes)+'\n\n这些散点按预先统一的覆盖及效应大小规则选取，为探索性展示。全组合(包括相反方向、接近零和不可估计)见该片 Pearson/delta 表。常规 Pearson P 不处理空间自相关；缺少空间区间不能解读为没有相关。正常来源不定义肿瘤 P95。\n';(O/'reports'/f'{sec}.md').write_text(text)
 return sec+' report complete'
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as ex:
  for s in ex.map(report,sys.argv[1:] or M.section_id.tolist()):print(s,flush=True)
