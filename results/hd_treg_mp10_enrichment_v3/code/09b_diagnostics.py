from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/home/data/t070721/.local/share/fonts/corefonts-arial/Arial.TTF');plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.linewidth':1.13,'lines.linewidth':1.3,'pdf.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
FD=R/'figures';idx=[]
def save(f,a,n,desc):
 f.savefig(FD/(n+'.png'),dpi=300,bbox_inches='tight');f.savefig(FD/(n+'.pdf'),bbox_inches='tight');leg=a.get_legend()
 if leg:
  h,l=a.get_legend_handles_labels();lf,la=plt.subplots(figsize=(5,1.3));la.axis('off');la.legend(h,l,loc='center',frameon=False);lf.savefig(FD/(n+'_legend_only.pdf'),bbox_inches='tight');plt.close(lf);leg.remove();f.savefig(FD/(n+'_panel_only.pdf'),bbox_inches='tight')
 plt.close(f);idx.append(dict(figure=n,description=desc,PNG=n+'.png',source='tables / compact regional data'))
def rd(f):return pd.read_csv(R/'tables'/f,sep='\t')
mi=rd('mixing_scores.tsv');hist=rd('null_histograms.tsv');ne=rd('neighborhood_enrichment.tsv');bo=rd('bootstrap_diagnostics.tsv');rob=rd('spatial_robustness.tsv');ef=rd('enrichment_effects.tsv');qq=[]
for s in SECTIONS:
 n=pd.read_parquet(R/s/'neighborhoods.parquet');u=n[n.radius_um.eq(200)&n.U_r];f,a=plt.subplots(figsize=(6.5,4.6));a.scatter(u.Treg_main,u.Epi_MP10,s=9,c='#6c4a8d',alpha=.5,rasterized=True);a.set_xlabel('Treg-related program (UCell)');a.set_ylabel('Independent Epi64 mean MP10');save(f,a,s+'_continuous_r200','All U200; no strict/FOXP3 gate. One T32 per point; overlapping neighborhoods not independent patients.')
 for radius in RADII:
  h=hist[hist.section.eq(s)&hist.radius_um.eq(radius)&hist.label_definition.eq('FOXP3_supported')&hist.high_definition.eq('H05')&hist.statistic.eq('edge_count')];z=ne[ne.section.eq(s)&ne.radius_um.eq(radius)&ne.label_definition.eq('FOXP3_supported')&ne.high_definition.eq('H05')].iloc[0];f,a=plt.subplots(figsize=(6.5,4.6));a.bar(h.null_value,h.frequency/9999,color='#c7b8d7',width=.85,label='Restricted-label null');a.axvline(z.O,color='#9c3985',label='Observed edges');a.set_xlabel('FOXP3-supported / H05 edges');a.set_ylabel('Null frequency');a.legend(frameon=False);save(f,a,s+'_null_edges_r'+str(radius),'9999 samples; counts retained even where no valid formal inference, e.g. zero null SD.')
 t=pd.read_parquet(R/s/'T_regions.parquet');t=t[t.T_eligible]
 for name,mask in [('FOXP3_supported',t.FOXP3_supported),('other_T',~t.FOXP3_supported),('strict',t.strict),('FOXP3_only',t.FOXP3_supported&~t.strict),('strict_only',t.strict&~t.FOXP3_supported)]:
  q=t[mask];qq.append(dict(section=s,group=name,n=len(q),zero_program=int(q.Treg_main.eq(0).sum()),median_program=q.Treg_main.median(),median_nFeature=q.nFeature.median(),median_nCount=q.nCount.median(),less3_signature_genes=int(q.Treg_main_detected.lt(3).sum()),no_core=int(q.Treg_core_detected.lt(1).sum()),below_old_cutoff=int(q.Treg_main.le(.19955389201641047).sum())))
for radius in RADII:
 for direction in ['label_to_high_Epi','high_Epi_to_label']:
  f,a=plt.subplots(figsize=(8,5));d=mi[mi.radius_um.eq(radius)&mi.high_definition.eq('H05')&mi.direction.eq(direction)]
  for j,label in enumerate(LABELS):
   q=d[d.label_definition.eq(label)].set_index('section').loc[SECTIONS];ys=np.arange(5)+(j-.5)*.2;a.scatter(q.MS,ys,color=['#9c3985','#187e8b'][j],label=label)
   for y,v in zip(ys,q.MS):
    if not np.isfinite(v):a.text(.01,y,'NA',transform=a.get_yaxis_transform(),fontsize=9)
  a.set_xscale('symlog',linthresh=.1);a.set_yticks(range(5),SECTIONS);a.invert_yaxis();a.set_xlabel('MS (symlog; region mixing)');a.legend(frameon=False);save(f,a,f'MS_{direction}_r{radius}','Unordered A-A; cross pair once, isolates included. No pseudocount.')
f,a=plt.subplots(figsize=(8,5));
for j,size in enumerate([400,800]):
 q=bo[bo.radius_um.eq(200)&bo.high_definition.eq('H05')&bo.label_definition.eq('FOXP3_supported')&bo.block_um.eq(size)].set_index('section').loc[SECTIONS]
 for i,(s,row)in enumerate(q.iterrows()):
  y=i+(j-.5)*.2;a.plot(row.FE,y,'o',color=['#9c3985','#267ca3'][j],label=str(size)+' μm blocks'if i==0 else None)
  if np.isfinite(row.FE_low):a.hlines(y,row.FE_low,row.FE_high,color=['#9c3985','#267ca3'][j])
  else:a.text(.01,y,'CI NA',transform=a.get_yaxis_transform(),fontsize=8)
a.axvline(1,color='gray',ls='--');a.set_yticks(range(5),SECTIONS);a.invert_yaxis();a.set_xlabel('FE (FOXP3-supported / H05, 200 μm)');a.legend(frameon=False);save(f,a,'FE_block_sensitivity','Matched center-based block intervals; invalid-draw fractions in table.')
f,a=plt.subplots(figsize=(8,5))
for i,s in enumerate(SECTIONS):
 d=rob[rob.section.eq(s)&rob.radius_um.eq(200)&rob.high_definition.eq('H05')&rob.label_definition.eq('FOXP3_supported')&rob.removed_centers.gt(0)];v=d.FE[np.isfinite(d.FE)];v0=ef[ef.section.eq(s)&ef.radius_um.eq(200)&ef.high_definition.eq('H05')&ef.label_definition.eq('FOXP3_supported')].FE.iloc[0]
 if len(v):a.hlines(i,v.min(),v.max(),color='#9c3985',linewidth=2)
 a.plot(v0,i,'o',color='black');a.text(1.02,i,str(len(v))+' deletions',transform=a.get_yaxis_transform(),va='center',fontsize=10)
a.axvline(1,color='gray',ls='--');a.set_yticks(range(5),SECTIONS);a.invert_yaxis();a.set_xlabel('FE range after leaving out each positive hotspot');save(f,a,'hotspot_influence','Dots original; line min/max from prespecified64µm linked-positive components. Unchanged remaining Epi neighborhoods; no model selection.')
put(pd.DataFrame(qq),'label_measurement_diagnostics.tsv');old=rd('figure_index.tsv');put(pd.concat([old[~old.figure.isin([i['figure']for i in idx])],pd.DataFrame(idx)],ignore_index=True),'figure_index.tsv');print('DIAGNOSTIC_PLOTS',len(idx))
