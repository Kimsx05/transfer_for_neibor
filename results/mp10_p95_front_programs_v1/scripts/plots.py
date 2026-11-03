from analyze import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from PIL import Image
import subprocess,json,gzip
for sty in ['Regular','Bold']:font_manager.fontManager.addfont(subprocess.check_output(['fc-match','-f','%{file}','Arial:style='+sty],text=True))
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.titlesize':12,'axes.labelsize':12,'pdf.fonttype':42,'lines.linewidth':1.2,'axes.linewidth':1.1,'figure.facecolor':'white','axes.facecolor':'white','axes.spines.right':False,'axes.spines.top':False})
PAL=pd.read_csv(O/'tables/niche_palette.tsv',sep='\t').set_index('accepted_label').color.to_dict()
FOCUS=['HALLMARK_'+x for x in ['E2F_TARGETS','G2M_CHECKPOINT','TGF_BETA_SIGNALING','EPITHELIAL_MESENCHYMAL_TRANSITION','HYPOXIA','ANGIOGENESIS','TNFA_SIGNALING_VIA_NFKB','IL6_JAK_STAT3_SIGNALING','INFLAMMATORY_RESPONSE','INTERFERON_ALPHA_RESPONSE','INTERFERON_GAMMA_RESPONSE','UNFOLDED_PROTEIN_RESPONSE','REACTIVE_OXYGEN_SPECIES_PATHWAY','DNA_REPAIR','MTORC1_SIGNALING','PI3K_AKT_MTOR_SIGNALING','ANDROGEN_RESPONSE']]+[PE,'FLORES_CRISPR_PRO_INVASION_10']
MAIN=sorted(G[G.panel=='main'].gs_name.unique());index=[]
def short(x):return x.replace('HALLMARK_','').replace('REACTOME_','Reactome ').replace('GSE15659_RESTING_VS_ACTIVATED_TREG_','Treg resting vs activated ').replace('EPITHELIAL_MESENCHYMAL_TRANSITION','EMT').replace('__NO_MP10',' (no MP10 genes)').replace('_',' ')
def export(fig,name,keyaxes=(),legendartists=(),kind='summary',section=''):
 p=O/'figures'/name;p.parent.mkdir(parents=True,exist_ok=True);fig.savefig(str(p)+'.pdf',facecolor='white',bbox_inches='tight');
 if os.environ.get('MP10_PDF_ONLY')!='1':fig.savefig(str(p)+'.png',dpi=300,facecolor='white',bbox_inches='tight')
 if keyaxes or legendartists:
  for ax in keyaxes:ax.set_visible(False)
  for a in legendartists:a.set_visible(False)
  fig.savefig(str(p)+'_panel_only.pdf',bbox_inches='tight',facecolor='white')
  for ax in fig.axes:
   if ax not in keyaxes:ax.set_visible(False)
  for ax in keyaxes:ax.set_visible(True)
  for a in legendartists:a.set_visible(True)
  for t in fig.texts:t.set_visible(False)
  # For legends attached to hidden axes, put the legend into a separate figure instead.
  if legendartists and not keyaxes:
   handles,labels=[],[]
   for a in legendartists:
    handles+=list(a.legendHandles);labels+=[t.get_text() for t in a.texts]
   lf=plt.figure(figsize=(8,max(1,len(labels)*.3)));lf.legend(handles,labels,loc='center',frameon=False,ncol=2);lf.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight',facecolor='white');plt.close(lf)
  else:fig.savefig(str(p)+'_legend_only.pdf',bbox_inches='tight',facecolor='white')
 plt.close(fig);index.append(dict(name=name,section_id=section,kind=kind,pdf=str(p.relative_to(O))+'.pdf',png=str(p.relative_to(O))+'.png'))
def heat(mat,name,label,limit=.3,cmap='RdBu_r'):
 mat=mat.reindex(index=[p for p in MAIN if p in mat.index]+[p for p in mat.index if p not in MAIN])
 f,a=plt.subplots(figsize=(max(7,len(mat.columns)*1.35+5),max(4,len(mat)*.24+1.7)))
 cm=plt.get_cmap(cmap).copy();cm.set_bad('#dddddd');im=a.pcolormesh(np.ma.masked_invalid(mat.to_numpy(float)),cmap=cm,vmin=-limit,vmax=limit,edgecolors='white',linewidth=.1)
 a.set_yticks(np.arange(len(mat))+.5,[short(p) for p in mat.index]);a.set_xticks(np.arange(len(mat.columns))+.5,mat.columns,rotation=45,ha='right');a.invert_yaxis();a.set_aspect('auto');cb=f.colorbar(im,ax=a,pad=.03);cb.set_label(label);cb.solids.set_rasterized(False);export(f,name,[cb.ax])
def section_plot(sec):
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');r=M.set_index('section_id').loc[sec];sc=json.loads((O/'tables/plot_scales.json').read_text());leg=[]
 f,axs=plt.subplots(2,4,figsize=(21,10.5));a=axs.ravel();f.suptitle(sec+' | frozen v3.4 / Epi MP10 and literature programs',fontsize=15)
 for ax in a:ax.set_aspect('equal');ax.set_xlim(d.x.min()-5,d.x.max()+5);ax.set_ylim(d.y.max()+5,d.y.min()-5);ax.axis('off')
 a[0].imshow(Image.open(r.image_path));a[0].set_title('H&E / observed outer boundary')
 gp=R/'results/mp10_pathology_regions_v3_4_final/by_dataset'/r.dataset_id/sec/'geometry.json.gz'
 if gp.exists():
  with gzip.open(gp,'rt') as h:geom=json.load(h)
  ed=pd.DataFrame(geom['roi_edges']['data'],columns=geom['roi_edges']['columns']);ed=ed[ed.boundary_type=='candidate_observed_outer'];a[0].add_collection(LineCollection(ed[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2),colors='#ff7f00',linewidths=1.1))
 for j,(col,title,mask,cmap,lim) in enumerate([(RAW,'Mixed-spot MP10 AUCell / P95',d['query'],'viridis',(0,sc[RAW])),(PE,'Puram p-EMT Top30 / active front',np.ones(len(d),bool),'viridis',(0,sc[PE])),('FLORES_CRISPR_PRO_INVASION_10','Flores positive 10',np.ones(len(d),bool),'viridis',(0,sc['FLORES_CRISPR_PRO_INVASION_10'])),('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','Hallmark EMT (mixed spot)',np.ones(len(d),bool),'viridis',(0,sc['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'])),('d_front','Strict v3.4 signed distance',np.ones(len(d),bool),'RdBu_r',(-1000,1000)),('d_MP10','Nearest MP10-P95 distance',np.ones(len(d),bool),'viridis',(0,1000))],start=1):
  ax=a[j];ax.scatter(d.x,d.y,s=3,c='#ddd',linewidths=0);ok=np.asarray(mask)&d[col].notna().to_numpy();im=ax.scatter(d.x[ok],d.y[ok],c=d.loc[ok,col],s=4,cmap=cmap,vmin=lim[0],vmax=lim[1],linewidths=0);cb=f.colorbar(im,ax=ax,shrink=.7,pad=.01);cb.set_label('μm' if col.startswith('d_') else 'AUCell' if col==RAW else 'UCell');cb.solids.set_rasterized(False);leg.append(cb.ax);ax.set_title(title)
  sel=d.MP10_P95 if col==RAW else d.pEMT_high_front if col==PE else np.zeros(len(d),bool)
  if np.any(sel):ax.scatter(d.x[sel],d.y[sel],s=12,facecolors='none',edgecolors='#D73027',linewidths=1)
 a[7].scatter(d.x,d.y,s=4,c=[PAL.get(v,'#ddd') for v in d.accepted_label.astype(object)],linewidths=0);a[7].set_title('SAMPLE R12 LEIDEN 0–8')
 lx=f.add_axes([.1,.005,.8,.035]);lx.axis('off');lx.legend(handles=[Patch(color=PAL[k],label='R12 '+str(k)) for k in range(9)]+[Line2D([],[],marker='o',color='#D73027',fillstyle='none',linestyle='',label='P95 / pEMT-high front')],loc='center',ncol=10,frameon=False,fontsize=10);leg.append(lx)
 f.subplots_adjust(left=.025,right=.98,top=.92,bottom=.08,wspace=.25,hspace=.28);name='sections/'+r.dataset_id+'/'+sec;export(f,name,leg,kind='section',section=sec)
 return index[-1]
def aggregate_curves():
 frames=[];rings=[]
 for s in M.section_id:
  frames.append(pd.read_csv(O/'tables'/f'{s}_curves.tsv.gz',sep='\t'));rings.append(pd.read_csv(O/'tables'/f'{s}_rings.tsv.gz',sep='\t'))
 def agg(a,keys):
  a=a.merge(M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid','eligible_for_patient_balanced_analysis']],on='section_id')
  outputs=[];dsets=[]
  for source,g in a.groupby('source_context'):
   for pop,z,idcol in [('sections',g,'section_id'),('units',g[g.analysis_unit_uid.notna()],'analysis_unit_uid'),('verified',g[g.eligible_for_patient_balanced_analysis.eq(True)&g.patient_uid.notna()],'patient_uid')]:
    val=z.dropna(subset=['mean']);u=val.groupby(keys+['dataset_id',idcol],dropna=False)[['mean','distance']].mean().reset_index();ds=u.groupby(keys+['dataset_id'],dropna=False)[['mean','distance']].mean().reset_index();ds['source']=source;ds['population_level']=pop;dsets.append(ds)
    sums=ds.groupby(keys,dropna=False).agg(mean=('mean','mean'),distance=('distance','mean'),n_datasets=('dataset_id','nunique'),se=('mean','sem')).reset_index()
    cnt=val.groupby(keys,dropna=False).agg(n_sections=('section_id','nunique'),n_units=(idcol,'nunique'),n_spots=('n','sum')).reset_index()
    vc=val[val.eligible_for_patient_balanced_analysis.eq(True)].groupby(keys,dropna=False).patient_uid.nunique().reset_index(name='n_verified_patients')
    sums=sums.merge(cnt,on=keys).merge(vc,on=keys,how='left');sums['source']=source;sums['population_level']=pop;outputs.append(sums)
  return pd.concat(outputs,ignore_index=True),pd.concat(dsets,ignore_index=True)
 c,cd=agg(pd.concat(frames,ignore_index=True),['axis','reference','population','bin','program']);t,td=agg(pd.concat(rings,ignore_index=True),['reference','population','ring','program']);save(c,O/'tables/curve_summary.tsv.gz');save(cd,O/'tables/curve_dataset.tsv.gz');save(t,O/'tables/ring_summary.tsv.gz');return c,t

def main():
 S=pd.read_csv(O/'tables/summary_effects.tsv',sep='\t');T=S[(S.source=='tumor_study_specimen_not_region_diagnosis')&(S.population=='units')];cov=pd.read_csv(O/'full_sections.tsv',sep='\t')
 c,t=aggregate_curves();c=c[(c.source=='tumor_study_specimen_not_region_diagnosis')&(c.population_level=='units')];t=t[(t.source=='tumor_study_specimen_not_region_diagnosis')&(t.population_level=='units')]
 for mode in ['raw','adjusted']:
  zt=T[(T.adjustment==mode)&T.analysis.isin(['B1_'+x for x in RINGS])&T.program.isin(MAIN)];heat(zt.pivot(index='program',columns='analysis',values='effect').reindex(columns=['B1_'+x for x in RINGS]),'02_center_rings_'+mode,'Standardized MP10 coefficient ('+mode+')',.4)
  zz=T[(T.adjustment==mode)&T.analysis.isin(['A_front_vs_internal','A_front_distance','C3_P95_state'])&T.program.isin(MAIN)];heat(zz.pivot(index='program',columns='analysis',values='effect'),'01_front_programs_'+mode,'Standardized effect ('+mode+')',.4)
 for typ in ['base','front']:
  zz=T[(T.adjustment=='matched')&T.analysis.str.startswith('B4_matched_'+typ+'_')&T.program.isin(MAIN)];heat(zz.pivot(index='program',columns='analysis',values='effect').reindex(columns=['B4_matched_'+typ+'_'+x for x in RINGS]),'02_matched_'+typ,'P95 minus matched control (section SD)',.5)
 zz=T[T.analysis.eq('B1_R0')&T.adjustment.eq('adjusted')];names=FOCUS+[p+'__NO_MP10' for p in FOCUS if p+'__NO_MP10' in P];zz=zz[zz.program.isin(names)].set_index('program').reindex(names).dropna(subset=['effect'])
 f,a=plt.subplots(figsize=(9,10));a.errorbar(zz.effect,np.arange(len(zz)),xerr=[zz.effect-zz.low,zz.high-zz.effect],fmt='o',color='#1f78b4',capsize=2);a.set_yticks(np.arange(len(zz)),[short(p) for p in zz.index]);a.invert_yaxis();a.axvline(0,color='grey',lw=1);a.set_xlabel('R0 adjusted MP10 coefficient / hierarchical 95% CI');export(f,'02_R0_overlap_sensitivity')
 # P95 field and paired-control distance curves, every prespecified non-lipid focus.
 for page in range(0,len(FOCUS),4):
  f,axs=plt.subplots(2,2,figsize=(14,10));leg=[]
  for ax,p in zip(axs.ravel(),FOCUS[page:page+4]):
   for pop,ls in [('tissue','-'),('Epi','--')]:
    for ref,col in [('P95_matched_front','#d73027'),('control_front','#1f78b4')]:
     zc=c[(c.axis=='d_MP10')&(c.program==p)&(c.population==pop)&(c.reference==ref)&c.bin.le(10)].set_index('bin').reindex(range(11)).reset_index();ax.plot(zc.distance,zc['mean'],ls,color=col,marker='o',ms=3,label=ref.replace('_',' ')+' / '+pop)
   ax.set_title(short(p));ax.set_xlabel('Distance to reference set (μm)');ax.set_ylabel('UCell / dataset-equal mean');ax.set_xlim(0,1000)
   zc=c[(c.axis=='d_MP10')&(c.program==p)&(c.population=='Epi')&(c.reference=='P95_matched_front')&c.bin.le(10)];ax.text(.02,.02,'Coverage: sections '+str(int(zc.n_sections.min()))+'–'+str(int(zc.n_sections.max()))+'; datasets '+str(int(zc.n_datasets.min()))+'–'+str(int(zc.n_datasets.max())) if len(zc) else 'No supported coverage',transform=ax.transAxes,fontsize=10)
  for ax in axs.ravel()[len(FOCUS[page:page+4]):]:ax.axis('off')
  handles,labels=axs.ravel()[0].get_legend_handles_labels();legax=f.add_axes([.1,.005,.8,.045]);legax.axis('off');legax.legend(handles,labels,loc='center',ncol=2,frameon=False,fontsize=10);f.subplots_adjust(hspace=.4,wspace=.3,bottom=.13);export(f,'03_P95_control_curves_'+str(page//4+1),[legax])
 # Front score curves in all tumor Epi versus P95, with outside tissue separate.
 for page in range(0,len([RAW]+FOCUS),4):
  f,axs=plt.subplots(2,2,figsize=(14,10))
  for ax,p in zip(axs.ravel(),([RAW]+FOCUS)[page:page+4]):
   for pop,col in [('tumor_Epi','#1f78b4'),('P95_Epi','#d73027'),('outside_tissue','#777777')]:
    zc=c[(c.axis=='d_front')&(c.program==p)&c.population.eq(pop)&c.bin.between(-1000,500)].set_index('bin').reindex(range(-1000,501,100)).reset_index();ax.plot(zc.distance,zc['mean'],color=col,marker='o',ms=3,label=pop);ax.fill_between(zc.distance,(zc['mean']-zc.se).to_numpy(),(zc['mean']+zc.se).to_numpy(),color=col,alpha=.1)
   ax.axvspan(-200,0,color='#ffe6a8',alpha=.4);ax.set_title(short(p));ax.set_xlabel('Strict v3.4 signed distance (μm)');ax.set_ylabel('AUCell' if p==RAW else 'UCell');zc=c[(c.axis=='d_front')&(c.program==p)&c.population.eq('tumor_Epi')&c.bin.between(-1000,0)];ax.text(.02,.02,'Sections per bin: '+str(int(zc.n_sections.min()))+'–'+str(int(zc.n_sections.max())) if len(zc) else 'No coverage',transform=ax.transAxes,fontsize=10)
  handles,labels=axs.ravel()[0].get_legend_handles_labels();la=f.add_axes([.1,.005,.8,.04]);la.axis('off');la.legend(handles,labels,loc='center',ncol=3,frameon=False);f.subplots_adjust(hspace=.4,wspace=.3,bottom=.13);export(f,'04_front_curves_'+str(page//4+1),[la])
 # Separate all-Epi/front MP10 effects from the P95 probability difference.
 f,axs=plt.subplots(1,2,figsize=(14,5))
 for ax,analyses,lab in [(axs[0],['C1_all_Epi_continuous','C1_front_vs_internal','C_active_front_MP10'],'Standardized MP10 coefficient / 95% CI'),(axs[1],['C2_P95_enrichment'],'P95 probability difference / 95% CI')]:
  zz=T[T.analysis.isin(analyses)&T.adjustment.eq('adjusted')].set_index('analysis').reindex(analyses).dropna(subset=['effect'])
  ax.errorbar(zz.effect,np.arange(len(zz)),xerr=[np.maximum(0,zz.effect-zz.low),np.maximum(0,zz.high-zz.effect)],fmt='o',capsize=3,color='#1f78b4');ax.set_yticks(np.arange(len(zz)),[x.replace('C1_','').replace('C2_','').replace('C_','').replace('_',' ') for x in zz.index]);ax.axvline(0,color='grey');ax.set_xlabel(lab)
 f.tight_layout();export(f,'04_MP10_front_effects')
 # Coverage is displayed independently, with sections / verified patients per bin.
 cc=c[(c.axis=='d_MP10')&c.program.eq(PE)&c.reference.isin(['P95_matched_front','control_front'])&c.bin.le(10)].copy();cc['row']=cc.reference+' / '+cc.population
 mat=cc.pivot(index='row',columns='bin',values='n_sections').reindex(columns=range(11));pat=cc.pivot(index='row',columns='bin',values='n_verified_patients').reindex(index=mat.index,columns=range(11))
 f,a=plt.subplots(figsize=(14,4));im=a.pcolormesh(np.ma.masked_invalid(mat.to_numpy(float)),cmap='Blues');a.set_xticks(np.arange(11)+.5,['0']+[str(i*100) for i in range(1,11)]);a.set_yticks(np.arange(len(mat))+.5,mat.index);a.set_xlabel('Reference-set distance bin upper bound (μm)')
 for i in range(len(mat)):
  for j in range(11):
   if np.isfinite(mat.iloc[i,j]):a.text(j+.5,i+.5,str(int(mat.iloc[i,j]))+'/'+str(int(pat.iloc[i,j])) if np.isfinite(pat.iloc[i,j]) else str(int(mat.iloc[i,j]))+'/0',ha='center',va='center',fontsize=10)
 cb=f.colorbar(im,ax=a);cb.set_label('Sections (cell text: sections / verified patients)');cb.solids.set_rasterized(False);export(f,'03_distance_coverage',[cb.ax])
 # Niche main effects all nine.
 zz=T[T.analysis.isin(['D_front_enrichment','D_active_front','D_P95_position','D_continuous_MP10_position','D_matched_front'])&T.adjustment.isin(['adjusted','matched'])];heat(zz.pivot(index='program',columns='analysis',values='effect'),'05_niche_effects','Niche probability difference / continuous MP10 SD',.15)
 # Coverage and P95 enrichment with correct denominator.
 b=pd.read_csv(O/'tables/bands.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','source_context','analysis_unit_uid']],on='section_id');b=b[(b.source_context!='adjacent_normal_source')&b.band.isin(['deep','near_inside','front'])];b['band']=pd.Categorical(b.band,['deep','near_inside','front'],ordered=True)
 f,a=plt.subplots(figsize=(10,5));order=['deep','near_inside','front']
 for ds,g in b.groupby('dataset_id'):
  v=g[g.analysis_unit_uid.notna()].groupby(['band','analysis_unit_uid'],observed=True).P95_fraction.mean().groupby('band',observed=True).mean().reindex(order);a.plot(range(3),v,marker='o',label=ds)
 a.set_xticks(range(3),order);a.set_ylabel('MP10-P95 / eligible tumor Epi');lg=a.legend(bbox_to_anchor=(1.01,1),frameon=False,fontsize=10);export(f,'04_P95_front_proportions',legendartists=[lg])
 n=pd.read_csv(O/'tables/niche.tsv.gz',sep='\t').merge(M[['section_id','dataset_id','analysis_unit_uid','source_context']],on='section_id');n=n[(n.source_context!='adjacent_normal_source')&n.analysis_unit_uid.notna()];order=['deep','near_inside','front','pEMT_high_front','other_front'];n=n[n.band.isin(order)]
 u=n.groupby(['dataset_id','analysis_unit_uid','band','niche']).fraction.mean().reset_index();ds=u.groupby(['dataset_id','band','niche']).fraction.mean().reset_index();v=ds.groupby(['band','niche']).fraction.mean().unstack().reindex(order);save(v.reset_index(),O/'tables/niche_band_equal_dataset.tsv')
 f,a=plt.subplots(figsize=(10,5));bottom=np.zeros(len(v))
 for k in range(9):a.bar(range(len(v)),v[k],bottom=bottom,color=PAL[k],label='R12 '+str(k));bottom+=v[k].fillna(0).to_numpy()
 a.set_xticks(range(len(v)),order,rotation=20,ha='right');a.set_ylabel('Fraction of labelled eligible Epi centers');lg=a.legend(bbox_to_anchor=(1.01,1),frameon=False);export(f,'05_niche_front_composition',legendartists=[lg])
 # Single-cell source summaries, sample medians equally weighted within type.
 scp=O/'tables/sc_sample_scores.tsv'
 if scp.exists():
  ss=pd.read_csv(scp,sep='\t');ss=ss.groupby(['sample','lineage','dataset'])[[PE,'FLORES_CRISPR_PRO_INVASION_10']].median().reset_index();v=ss.groupby('lineage')[[PE,'FLORES_CRISPR_PRO_INVASION_10']].median();f,axs=plt.subplots(1,2,figsize=(11,4));
  for ax,p in zip(axs,v.columns):
   groups=[ss.loc[ss.lineage.eq(x),p].dropna() for x in v.index];ax.boxplot(groups,labels=v.index,showfliers=False);ax.set_ylabel('Sample × subtype median UCell');ax.set_title(short(p));ax.tick_params(axis='x',rotation=25)
  f.tight_layout();export(f,'01_single_cell_source')
 # Predetermined shared scales from all observed values; display only, no effect on statistics.
 cols=[RAW,PE,'FLORES_CRISPR_PRO_INVASION_10','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'];values={k:[] for k in cols}
 for s in M.section_id:
  dd=pd.read_parquet(O/'spots'/f'{s}.parquet',columns=cols)
  for k in cols:values[k].extend(dd[k].dropna().tolist())
 sc={k:float(np.quantile(v,.98)) for k,v in values.items()};(O/'tables/plot_scales.json').write_text(json.dumps(sc,indent=2))
 summary_index=index.copy()
 with ProcessPoolExecutor(max_workers=4) as ex:
  for a in ex.map(section_plot,M.section_id):summary_index.append(a);print(a['section_id'],flush=True)
 save(summary_index,O/'figure_index.tsv')
if __name__=='__main__':main()
