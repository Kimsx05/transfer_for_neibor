from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors,cm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF')
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1.13,'lines.linewidth':1.3,'pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
FD=R/'figures';FD.mkdir(exist_ok=True);idx=[];C={'FOXP3_supported':'#9c3985','strict':'#187e8b','H05':'#d95f02','other':'#aeb8bf'}
def rd(n):return pd.read_csv(R/'tables'/n,sep='\t')
def save(fig,ax,name,desc,cb=None):
 
 only=os.environ.get('PLOT_UPDATE_NAMES','')
 if only and name not in only.split(','):
  plt.close(fig);return
 leg=ax.get_legend();handles,labels=ax.get_legend_handles_labels()
 if leg is not None:handles=leg.legend_handles;labels=[z.get_text()for z in leg.get_texts()]
 fig.savefig(FD/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(FD/(name+'.pdf'),bbox_inches='tight')
 if leg is not None or cb is not None:
  lf,la=plt.subplots(figsize=(5, max(1.2,.32*len(labels))));la.axis('off')
  if handles:la.legend(handles,labels,loc='center',frameon=False)
  if cb is not None:
   cf,ca=plt.subplots(figsize=(4.8,1.1));cf.subplots_adjust(bottom=.55,top=.8);cf.colorbar(cm.ScalarMappable(norm=cb.norm,cmap=cb.cmap),cax=ca,orientation='horizontal',label=cb.ax.get_ylabel());cf.savefig(FD/(name+'_legend_only.pdf'),bbox_inches='tight');plt.close(cf)
  else:lf.savefig(FD/(name+'_legend_only.pdf'),bbox_inches='tight')
  plt.close(lf)
  if leg is not None:leg.remove()
  if cb is not None:cb.remove()
  fig.savefig(FD/(name+'_panel_only.pdf'),bbox_inches='tight')
 plt.close(fig);idx.append(dict(figure=name,description=desc,PNG=name+'.png',source='tables and section compact Parquet; original members referenced read-only'))
def spatial(ax):ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (μm)');ax.set_ylabel('y (μm)')
def basefig(w=6.5,h=4.6):return plt.subplots(figsize=(w,h))
counts=rd('label_counts_and_overlap.tsv');ef=rd('enrichment_effects.tsv');ne=rd('neighborhood_enrichment.tsv');co=rd('continuous_associations.tsv');cur=rd('cross_distance_results.tsv');mix=rd('mixing_scores.tsv');boot=rd('bootstrap_diagnostics.tsv')
f,a=basefig(8,4.8);x=np.arange(5)
for j,(col,color,lab)in enumerate([('T_eligible','#b4bec6','All evaluable T regions'),('FOXP3_supported',C['FOXP3_supported'],'FOXP3-supported T regions'),('strict',C['strict'],'Historical strict Treg-like')]):a.bar(x+(j-1)*.25,counts[col],width=.24,color=color,label=lab);[a.text(xx+(j-1)*.25,yy+2,str(yy),ha='center',va='bottom',fontsize=9)for xx,yy in zip(x,counts[col])]
a.set_yscale('symlog',linthresh=3);a.set_xticks(x,SECTIONS,rotation=20);a.set_ylabel('Regions (symlog scale)');a.legend(frameon=False);save(f,a,'definition_counts','All five sections; T32 regions, not cells.')
f,a=basefig(8,4.8);bottom=np.zeros(5)
for col,color,lab in [('FOXP3_only',C['FOXP3_supported'],'FOXP3 only'),('intersection','#74449b','Both labels'),('strict_only',C['strict'],'Strict only')]:a.bar(x,counts[col],bottom=bottom,label=lab,color=color);bottom+=counts[col]
a.set_xticks(x,SECTIONS,rotation=20);a.set_ylabel('Label-positive regions');a.legend(frameon=False);save(f,a,'label_overlap','Mutually exclusive overlap categories; remaining eligible T regions in table.')
# Forest plots retain every planned row, including undefined effect or missing intervals.
for radius in RADII:
 f,a=basefig(8,5.3)
 for j,label in enumerate(LABELS):
  d=ef[(ef.radius_um==radius)&ef.high_definition.eq('H05')&ef.label_definition.eq(label)].set_index('section').loc[SECTIONS]
  for i,(s,row)in enumerate(d.iterrows()):
   y=i+(j-.5)*.22;v=row.FE
   if np.isfinite(v):
    a.plot(v,y,'o'if j==0 else's',color=C[label],label=label if i==0 else None)
    if np.isfinite(row.FE_low):a.hlines(y,row.FE_low,row.FE_high,color=C[label],linewidth=1.5)
    else:a.annotate('CI NA',(v,y),xytext=(5,-9),textcoords='offset points',fontsize=8,color=C[label])
   else:a.text(.02,y,'NA',transform=a.get_yaxis_transform(),color=C[label],fontsize=9)
   a.text(1.02,y,f'{int(row.label_positive)}/{int(row.N)}',transform=a.get_yaxis_transform(),va='center',fontsize=9,color=C[label])
 a.axvline(1,color='#777777',ls='--');a.set_yticks(x,SECTIONS);a.invert_yaxis();a.set_xlabel('FE (H05; unique T regions)');a.set_ylim(4.6,-.6);a.legend(handles=[Line2D([0],[0],color=C[l],marker='o'if l==LABELS[0]else's',label=l)for l in LABELS],frameon=False,loc='upper center',bbox_to_anchor=(.5,1.17),ncol=2);save(f,a,f'FE_H05_r{radius}','400μm spatial-block intervals; NA retained. Right: positive/all U_r centers. Strict is historical sensitivity.')
 f,a=basefig(8,5.3)
 for j,label in enumerate(LABELS):
  d=ne[(ne.radius_um==radius)&ne.high_definition.eq('H05')&ne.label_definition.eq(label)].set_index('section').loc[SECTIONS]
  for i,(s,row)in enumerate(d.iterrows()):
   y=i+(j-.5)*.22
   if np.isfinite(row.Z):a.plot(row.Z,y,'o',color=C[label],label=label if i==0 else None)
   else:a.text(.02,y,'NA',transform=a.get_yaxis_transform(),color=C[label],fontsize=9)
   a.text(1.02,y,f'O={int(row.O)} / T+={int(row.n_label)}',transform=a.get_yaxis_transform(),fontsize=9,va='center',color=C[label])
 a.axvline(0,color='#777777',ls='--');a.set_yticks(x,SECTIONS);a.invert_yaxis();a.set_xlabel('Edge enrichment Z (restricted label null)');a.legend(handles=[Line2D([0],[0],color=C[l],marker='o',label=l)for l in LABELS],frameon=False);save(f,a,f'edge_Z_H05_r{radius}','Edges weight centers by their number of high Epi neighbors; distinct from center-based FE.')
 for direction in ['label_to_high_Epi','high_Epi_to_label']:
  f,a=basefig(8,5.3)
  for j,label in enumerate(LABELS):
   d=mix[(mix.radius_um==radius)&mix.high_definition.eq('H05')&mix.label_definition.eq(label)&mix.direction.eq(direction)].set_index('section').loc[SECTIONS]
   for i,(s,row)in enumerate(d.iterrows()):
    y=i+(j-.5)*.22
    if np.isfinite(row.NMS):a.plot(row.NMS,y,'o',color=C[label],label=label if i==0 else None)
    else:a.text(.02,y,'NA',transform=a.get_yaxis_transform(),color=C[label],fontsize=9)
    a.text(1.02,y,f'{int(row.n_A)}/{int(row.n_B)}',transform=a.get_yaxis_transform(),va='center',fontsize=9,color=C[label])
  a.set_xscale('symlog',linthresh=.2);a.set_yticks(x,SECTIONS);a.invert_yaxis();a.set_xlabel('NMS (symlog scale; descriptive)');a.legend(handles=[Line2D([0],[0],color=C[l],marker='o',label=l)for l in LABELS],frameon=False);save(f,a,f'NMS_{direction}_r{radius}','Direction fixed in filename; right n_A/n_B includes isolated points. NMS=1 is not a tissue-adjusted null.')
 f,a=basefig(8,5.3)
 for j,size in enumerate([400,800]):
  d=co[(co.radius_um==radius)&co.readout.eq('Epi_MP10')&co.block_um.eq(size)].set_index('section').loc[SECTIONS]
  for i,(s,row)in enumerate(d.iterrows()):
   y=i+(j-.5)*.22;color=['#6a4a9a','#298e9a'][j];a.plot(row.adjusted_rank_beta,y,'o',color=color,label=f'{size} μm blocks'if i==0 else None)
   if np.isfinite(row.low):a.hlines(y,row.low,row.high,color=color)
   else:a.annotate('CI NA',(row.adjusted_rank_beta,y),xytext=(5,-8),textcoords='offset points',fontsize=8)
 a.axvline(0,color='gray',ls='--');a.set_yticks(x,SECTIONS);a.invert_yaxis();a.set_xlabel('Adjusted standardized rank effect');a.legend(frameon=False);save(f,a,f'continuous_effect_r{radius}','Continuous all eligible U_r, independent Epi mean; no new formal P. Prespecified covariates and 400/800μm blocks.')
# Calibration performance, denominators unchanged between old/fixed rules.
cal=rd('calibration_fixed_performance.tsv');ca=cal.groupby(['partition','rule','metric'],as_index=False)[['numerator','T_eligible_denominator']].sum();ca['rate']=ca.numerator/ca.T_eligible_denominator;put(ca,'calibration_performance_pooled.tsv')
f,a=basefig(8,4.8)
for j,rule in enumerate(['historical','denominator_fixed']):
 d=ca[ca.rule.eq(rule)];a.bar(np.arange(len(d))+(j-.5)*.3,d.rate,width=.29,label=rule,color=['#76828c','#9b498b'][j]);a.set_xticks(np.arange(len(d)),d.partition+'\n'+d.metric,rotation=15)
a.set_ylabel('Full-rule fraction among T-eligible mixtures');a.legend(frameon=False);save(f,a,'calibration_full_rule','Pooled display only; selection used every nonempty training condition. Holdout previously viewed. No HD labels replaced.')
roi=rd('roi_manifest.tsv');rois=[];plotlong=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet');u=t[t.T_eligible].copy();e=pd.read_parquet(R/s/'Epi_regions.parquet');n=pd.read_parquet(R/s/'neighborhoods.parquet');epall=pd.read_parquet(V/s/'Epi64_regions.parquet')
 for _,v in u.iterrows():plotlong.append(dict(section=s,group='FOXP3-supported'if v.FOXP3_supported else'Other evaluable T',Treg_main=v.Treg_main,log10_UMI=np.log10(v.nCount+1),log1p_FOXP3=np.log1p(v.raw_FOXP3),log1p_IL2RA=np.log1p(v.raw_IL2RA),log1p_CTLA4=np.log1p(v.raw_CTLA4)))
 # Full slice plots share exact extent within section. Background is supplied T32, not confirmed tissue.
 extent=[t.x_um.min()-32,t.x_um.max()+32,t.y_um.min()-32,t.y_um.max()+32]
 def background(a):a.scatter(t.x_um,t.y_um,s=1,c='#e4e7e9',rasterized=True);a.set_xlim(extent[:2]);a.set_ylim(extent[2:]);spatial(a)
 f,a=basefig(7,6);background(a);p=a.scatter(u.x_um,u.y_um,c=u.Treg_main,s=5,cmap='viridis',vmin=0,vmax=.55,rasterized=True);cb=f.colorbar(p,ax=a,label='Treg-related UCell (above0.55 clipped)',extend='max');save(f,a,s+'_program','All T-eligible scores including zero. Gray supplied T32 below measurement eligibility; white no represented input. MPBC relative coordinates, original tissue unknown.',cb)
 f,a=basefig(7,6);background(a)
 for mask,col,mark,lab in [(~u.FOXP3_supported&~u.strict,'#9caab3','.', 'Other evaluable T regions'),(u.FOXP3_supported&~u.strict,C['FOXP3_supported'],'o','FOXP3 only'),(~u.FOXP3_supported&u.strict,C['strict'],'s','Strict only'),(u.FOXP3_supported&u.strict,'#422667','D','Both labels')]:
  q=u[mask];a.scatter(q.x_um,q.y_um,s=8 if mark!='.'else 3,c=col,marker=mark,label=lab,rasterized=True)
 a.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(1,1));save(f,a,s+'_labels','Original strict unchanged; simple FOXP3 label has no UCell/core gate; background supplied ineligible T32.')
 f,a=basefig(7,6);background(a);a.scatter(epall.x_um,epall.y_um,s=3,c='#b8c4cb',rasterized=True);p=a.scatter(e.x_um,e.y_um,c=e.MP10,s=8,cmap='magma',vmin=0,vmax=.3,rasterized=True);cb=f.colorbar(p,ax=a,label='Qualified Epi64 MP10 (above0.3 clipped)',extend='max');save(f,a,s+'_Epi_MP10','Unqualified Epi regions gray-blue, not MP10-low; zero is a computed value; no score interpolation.',cb)
 f,a=basefig(7,6);background(a);a.scatter(e.x_um,e.y_um,s=4,c='#9caab3',label='Qualified Epi64, not H05',rasterized=True);q=e[e.H05];a.scatter(q.x_um,q.y_um,s=13,c=C['H05'],label='Epi64 H05',rasterized=True);a.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(1,1));save(f,a,s+'_Epi_H05','Within-section q95 strict greater-than; includes no tied-value forcing.')
 f,a=basefig(7,6);background(a);wide=n.pivot(index='region_id',columns='radius_um',values='U_r');uu=u.set_index('region_id').loc[wide.index]
 for mask,col,lab in [(~wide[50]&~wide[200],'#8c9196','Neither scale evaluable'),(~wide[50]&wide[200],'#4d91ad','200 μm only'),(wide[50]&~wide[200],'#e2b650','50 μm only'),(wide[50]&wide[200],'#743f8d','Both scales')]:a.scatter(uu.loc[mask,'x_um'],uu.loc[mask,'y_um'],s=7,c=col,label=lab,rasterized=True)
 a.legend(frameon=False,fontsize=9,loc='upper left',bbox_to_anchor=(1,1));save(f,a,s+'_evaluability','Full T-eligible denominator. Not evaluable is missing neighborhood measurement, not non-high exposure.')
 for label in LABELS:
  q=cur[cur.section.eq(s)&cur.high_definition.eq('H05')&cur.label_definition.eq(label)];f,a=basefig();a.plot(q.distance_um,q.observed_over_null,c=C[label],label='Observed / restricted-null mean');a.fill_between(q.distance_um,q.null_low_ratio,q.null_high_ratio,color=C[label],alpha=.2,label='Pointwise null envelope');a.axhline(1,color='gray',ls='--');a.axvline(50,color='#777777',ls=':');a.axvline(200,color='#777777',ls=':');a.set_xlim(0,200);a.set_xlabel('Cumulative centroid radius (μm)');a.set_ylabel('Source-excluded cross proximity ratio');a.legend(frameon=False,fontsize=9);save(f,a,s+'_distance_'+label,'Fixed U200 for every radius; same source mask and9999 restricted permutations; not Kcross or simultaneous envelope. Null support='+str(q.null_eligible.iloc[0]))
 # Fixed previous ROI centers: one readable overlay each, plus links to original raw marker atlas.
 pairs=pd.read_parquet(V/s/'region_pairs.parquet');pairs=pairs[pairs.T_grid_um.eq(32)&pairs.radius_um.eq(200)]
 for _,rr in roi[roi.section.eq(s)].iterrows():
  cx=(rr.x0+rr.x1)/2;cy=(rr.y0+rr.y1)/2;tt=t[(t.x_um-cx).abs()<=256];tt=tt[(tt.y_um-cy).abs()<=256];ee=e[(e.x_um-cx).abs()<=256];ee=ee[(ee.y_um-cy).abs()<=256];f,a=basefig(7,6)
  a.scatter(tt.x_um,tt.y_um,c='#e1e4e7',s=13,marker='s',label='Supplied T32 footprint centers');a.scatter(ee.x_um,ee.y_um,c=np.where(ee.H05,C['H05'],'#acb7bf'),s=26,marker='s',label='Epi64 (orange: H05)')
  for label in LABELS:
   zz=tt[tt[label]];a.scatter(zz.x_um,zz.y_um,s=85 if label==LABELS[0]else 32,marker='o'if label==LABELS[0]else'D',facecolors='none',edgecolors=C[label],linewidths=1.2,label=label)
  a.scatter([cx],[cy],marker='+',s=100,c='black',label='Fixed ROI T32 center')
  for r in RADII:a.add_patch(Circle((cx,cy),r,fill=False,ls='--'if r==50 else'-',edgecolor='black',linewidth=1.15));a.text(cx+r/1.414,cy-r/1.414,str(r)+' μm',fontsize=10)
  a.set_xlim(cx-256,cx+256);a.set_ylim(cy-256,cy+256);spatial(a);a.legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(1,1));name=s+'_ROI_'+rr.roi;save(f,a,name,'Same fixed ROI center;512μm display,50/200µm cumulative circles. Shows all nearby scoreable Epi incl pairs excluded from inference; pair counts in ROI table distinguish sources. No H&E fabricated; original raw evidence atlas referenced.')
  tr=t.set_index('region_id').loc[rr.region_id];out=dict(section=s,roi=rr.roi,region_id=rr.region_id,x_um=cx,y_um=cy,FOXP3_supported=bool(tr.FOXP3_supported),strict=bool(tr.strict),FOXP3_UMI=tr.raw_FOXP3,IL2RA_UMI=tr.raw_IL2RA,CTLA4_UMI=tr.raw_CTLA4,IKZF2_UMI=tr.raw_IKZF2,Treg_program=tr.Treg_main,Tfh=tr.Tfh,Activated_T=tr.Activated_T)
  for r in RADII:
   z=n[n.region_id.eq(rr.region_id)&n.radius_um.eq(r)];p=pairs[pairs.T_region_id.eq(rr.region_id)&pairs.centroid_distance_um.le(r)]
   if len(z):out.update({f'r{r}_U':bool(z.U_r.iloc[0]),f'r{r}_Epi_n':z.Epi_n.iloc[0],f'r{r}_H05_n':z.n_H05.iloc[0],f'r{r}_shared_pairs_removed':int(p.shared_original_bins.sum()),f'r{r}_MP10mean':z.Epi_MP10.iloc[0]})
  rois.append(out)
 print('PLOTS',s,flush=True)
# Vector descriptive K/L only where estimated; all NA component rows remain in table.
k=rd('standard_cross_KL.tsv')
for s in SECTIONS[:3]:
 q=k[k.section.eq(s)&k.high_definition.eq('H05')&k.label_definition.eq('FOXP3_supported')];f,a=basefig()
 for comp,z in q.groupby('component'):
  if z.L_minus_r.notna().any():a.plot(z.distance_um,z.L_minus_r,label='Input component '+str(comp))
 a.axhline(0,color='gray',ls='--');a.axvline(50,color='gray',ls=':');a.axvline(200,color='gray',ls=':');a.set_xlabel('Radius (μm)');a.set_ylabel('Border-estimate L(r) - r (μm)');a.legend(frameon=False,fontsize=9);save(f,a,s+'_Kcross_descriptive','Globally source-disjoint subset;8μm supplied in-tissue pixel window; homogeneous border estimate. Nonuniformity unresolved, not formal attraction test.')
put(pd.DataFrame(rois),'fixed_ROI_v3_evidence.tsv');pd.DataFrame(plotlong).to_parquet(R/'tables/distribution_plot_data.parquet',index=False);pd.DataFrame(plotlong).to_csv(R/'tables/distribution_plot_data.tsv.gz',sep='\t',index=False,compression='gzip');old=pd.read_csv(R/'tables/figure_index.tsv',sep='\t')if (R/'tables/figure_index.tsv').exists()else pd.DataFrame();new=pd.DataFrame(idx);put(pd.concat([old[~old.figure.isin(new.figure)]if len(old)else old,new],ignore_index=True),'figure_index.tsv')
print('ALL_PLOTS',len(idx),flush=True)
