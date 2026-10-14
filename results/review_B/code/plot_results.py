from common import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import textwrap,argparse
plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'pdf.fonttype':42,'savefig.transparent':True})
LEGEND=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#2166ac',markeredgecolor='none',label='Verified patient unit',markersize=4),Line2D([0],[0],marker='^',color='none',markerfacecolor='#d88924',markeredgecolor='none',label='Other policy unit',markersize=4)]
def boxes(ax,data,value,scale=1):
 vals=[]
 for ri,r in enumerate(RINGS):
  g=data[data.ring.eq(r)];x=g[value].dropna().to_numpy()*scale;vals.append(x)
  if len(x):
   ax.boxplot([x],positions=[ri],widths=.48,whis=1.5,showfliers=False,patch_artist=True,boxprops=dict(facecolor=['#dde6ef','#e5e1ec','#f3e3d5'][ri],edgecolor='#404040'),medianprops=dict(color='black',linewidth=1),whiskerprops=dict(color='#555555'),capprops=dict(color='#555555'))
   for known,marker,color,shift in [(True,'o','#2166ac',-.06),(False,'^','#d88924',.06)]:
    z=g[g.known_patient.eq(known)][value].to_numpy()*scale;ax.scatter(np.full(len(z),ri+shift),z,s=8,alpha=.5,marker=marker,c=color,edgecolors='none',zorder=3)
  else:ax.text(ri,.5,'N/A',transform=ax.get_xaxis_transform(),ha='center',color='#666666')
 ax.set_xticks(range(3),[f'{r}\nN={len(x)}' for r,x in zip(RINGS,vals)]);ax.set_xlim(-.55,2.55);ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='y',labelsize=7)
 if not any(len(x) for x in vals):ax.set_ylim(-.1,.1)
 return vals

def star_above(ax,qs):
 lo,hi=ax.get_ylim();span=max(hi-lo,1e-5);ax.set_ylim(lo,hi+span*.24)
 for ri,q in enumerate(qs):ax.text(ri,hi+span*.07,stars(q),ha='center',fontsize=8)

def brackets(ax,qs):
 lo,hi=ax.get_ylim();span=max(hi-lo,1e-5);ax.set_ylim(lo,hi+span*.65)
 for j,((a,b),q) in enumerate(zip([(0,1),(0,2),(1,2)],qs)):
  y=hi+span*(.08+.19*j);ax.plot([a,a,b,b],[y,y+.04*span,y+.04*span,y],lw=.65,color='#444444');ax.text((a+b)/2,y+.05*span,stars(q),ha='center',va='bottom',fontsize=7)

def save(fig,path):
 path.parent.mkdir(parents=True,exist_ok=True);fig.savefig(path.with_suffix('.png'),dpi=300,transparent=True);fig.savefig(path.with_suffix('.pdf'),transparent=True);plt.close(fig)

def plot_scope(scope):
 root=OUT/'06_figures'/scope;root.mkdir(parents=True,exist_ok=True);pp=pairs();order={x:i for i,x in enumerate(TIERS)};pp['order']=pp.review_tier.map(order);pp=pp.sort_values(['order','accepted_label','cell_index']);mapping=[]
 comp=pd.read_csv(OUT/'02_background_composition/unit_ring123_composition.tsv.gz',sep='\t');comp=comp[comp.scope.eq(scope)&comp.analysis_set.eq('MAIN')&comp.summary_type.eq('INFERENTIAL')];cz=pd.read_csv(OUT/'03_ring_statistics/ring123_background_tests.tsv',sep='\t');cz=cz[cz.scope.eq(scope)&cz.unit_subset.eq('ALL_POLICY_UNITS')];cp=pd.read_csv(OUT/'03_ring_statistics/ring123_pairwise_tests.tsv',sep='\t');cp=cp[cp.scope.eq(scope)&cp.unit_subset.eq('ALL_POLICY_UNITS')]
 cg={k:g for k,g in comp.groupby('pair_id')};zg={k:g.set_index('ring') for k,g in cz.groupby('pair_id')};pg={k:g.set_index('ring_pair') for k,g in cp.groupby('pair_id')}
 for n,pr in enumerate(pp.itertuples()):
  rel='composition/'+pr.pair_id+'_excess';path=root/rel;g=cg.get(pr.pair_id,comp.iloc[:0]);z=zg[pr.pair_id];pair=pg[pr.pair_id]
  if not path.with_suffix('.png').exists() or not path.with_suffix('.pdf').exists():
   fig,axes=plt.subplots(1,3,figsize=(10.5,4.3));fig.subplots_adjust(left=.07,right=.98,bottom=.28,top=.72,wspace=.45)
   boxes(axes[0],g,'target',100);axes[0].set_title('Raw target composition\nNo inferential stars');axes[0].set_ylabel('Composition (%)')
   boxes(axes[1],g,'excess',100);axes[1].axhline(0,c='#777777',ls='--',lw=.7);star_above(axes[1],[z.loc[r,'background_zero_q'] for r in RINGS]);axes[1].set_title('Target minus same-section background\nq: background_zero_q');axes[1].set_ylabel('Difference (percentage points)')
   boxes(axes[2],g,'excess',100);axes[2].axhline(0,c='#777777',ls='--',lw=.7);brackets(axes[2],[pair.loc[r,'ring_pairwise_q'] for r in ['R1-R2','R1-R3','R2-R3']]);axes[2].set_title('Paired ring differences\nq: ring_pairwise_q');axes[2].set_ylabel('Difference (percentage points)')
   fig.suptitle(scope+' | '+pr.niche+' | '+pr.celltype+'\n'+pr.review_tier+' (selection tier; does not determine test direction)',fontsize=10,y=.97)
   fig.legend(handles=LEGEND,loc='lower center',bbox_to_anchor=(.5,.09),ncol=2,frameon=False,fontsize=8)
   fig.text(.07,.025,'Two-sided unit signed-rank; independent BH families of 417 per scope. N/A = insufficient units/nonzero pairs.\nSame eligible sections and units in all rings; boxes: median, Q1/Q3, 1.5 IQR; all unit points shown once. Exploratory.',fontsize=7)
   save(fig,path)
  mapping.append(dict(scope=scope,kind='composition',query_id=pr.pair_id,niche=pr.niche,celltype=pr.celltype,review_tier=pr.review_tier,png=f'06_figures/{scope}/{rel}.png',pdf=f'06_figures/{scope}/{rel}.pdf',n_units=int(z.n_units.max())))
  if n%35==0:log(scope+f' composition {n+1}/139')
 # One association figure has 2 score layers x main/raw panels; no duplicate ALL query per niche.
 au=pd.read_csv(OUT/'04_MP10_association/MP10_unit_associations.tsv.gz',sep='\t');au=au[au.scope.eq(scope)&au.analysis_set.eq('MAIN')&au.method.isin(['ADJUSTED_MAIN','RAW'])];at=pd.read_csv(OUT/'04_MP10_association/MP10_scope_associations.tsv',sep='\t');at=at[at.scope.eq(scope)&at.unit_subset.eq('ALL_POLICY_UNITS')];ag={k:g for k,g in au.groupby(['query_id','method','score_layer'])};tg={k:g.set_index('ring') for k,g in at.groupby(['query_id','method','score_layer'])}
 for n,query in enumerate(queries()):
  rel='MP10/'+query['query_id']+'_association';path=root/rel
  if not path.with_suffix('.png').exists() or not path.with_suffix('.pdf').exists():
   fig,axes=plt.subplots(2,2,figsize=(8.3,6.8));fig.subplots_adjust(left=.10,right=.97,bottom=.18,top=.79,hspace=.72,wspace=.32)
   for row,layer in enumerate(LAYERS):
    for col,method in enumerate(['ADJUSTED_MAIN','RAW']):
     ax=axes[row,col];key=(query['query_id'],method,layer);g=ag.get(key,au.iloc[:0]);t=tg[key];boxes(ax,g,'rho');ax.axhline(0,c='#777777',ls='--',lw=.7);ax.set_ylim(-1,1);star_above(ax,[t.loc[r,'association_q'] for r in RINGS]);ax.set_ylabel('Unit correlation');ax.set_title(('Adjusted main' if method=='ADJUSTED_MAIN' else 'Raw Spearman (auxiliary)')+' | '+('whole spot' if row==0 else 'model-attributed Epi')+'\nq: association_q / '+method)
   label='all eligible epithelial centers' if query['accepted_label']==-99 else f"within niche{query['accepted_label']}"
   fig.suptitle(scope+' | '+label+'\n'+query['celltype']+' | continuous MP10 association',fontsize=10,y=.97);fig.legend(handles=LEGEND,loc='lower center',bbox_to_anchor=(.5,.085),ncol=2,frameon=False,fontsize=8)
   fig.text(.10,.025,'Two-sided unit signed-rank vs 0; BH planned n=1008 per method/scope. Main and raw q are separate.\nSection effects -> equal-section Fisher-z average per unit. Same scored centers and paired score-layer sections.\nN/A: insufficient units/nonzero pairs. No spot-level P used; model-attributed Epi is not measured pure Epi.',fontsize=7)
   save(fig,path)
  mapping.append(dict(scope=scope,kind='MP10',query_id=query['query_id'],niche='ALL_ELIGIBLE_EPI' if query['accepted_label']==-99 else f"niche{query['accepted_label']}",celltype=query['celltype'],review_tier='SHARED_ALL_RANGE' if query['accepted_label']==-99 else pp.set_index('pair_id').loc[query['query_id'],'review_tier'],png=f'06_figures/{scope}/{rel}.png',pdf=f'06_figures/{scope}/{rel}.pdf',n_units=int(at[at.query_id.eq(query['query_id'])].n_units.max())))
  if n%40==0:log(scope+f' MP10 {n+1}/168')
 tab(f'06_figures/{scope}/display_name_map.tsv',mapping)
 # Structural Epi/TME summaries: nine niches, three adjacent boxes each, no new mechanistic tests.
 cl=pd.read_csv(OUT/'05_sensitivity/closure_unit_summary.tsv',sep='\t');cl=cl[cl.scope.eq(scope)]
 fig,axes=plt.subplots(3,3,figsize=(10,9));fig.subplots_adjust(top=.88,bottom=.1,hspace=.6,wspace=.35)
 for k,ax in enumerate(axes.flat):
  g=cl[cl.niche.eq(f'niche{k}')];boxes(ax,g,'Epi_ring',100);ax.set_ylim(0,100);ax.set_ylabel('Epithelial share (%)');ax.set_title(f'niche{k}: TME = 100% - Epi')
 fig.suptitle(scope+' | structural composition (descriptive units)\nniche0 remains in background; no niche0 celltype mechanism branches',fontsize=11);fig.legend(handles=LEGEND,loc='lower center',ncol=2,frameon=False);save(fig,root/'structural_Epi_TME')
 # Compact overview uses numeric effect matrices; annotations do not replace full boxes.
 fig,axes=plt.subplots(1,3,figsize=(12,9));fig.subplots_adjust(left=.27,right=.94,bottom=.08,top=.9,wspace=.4)
 core=pp[pp.review_tier.eq('POSITIVE_CORE')]
 for ax,ring in zip(axes,RINGS):
  Z=cz[cz.ring.eq(ring)].set_index('pair_id').reindex(core.pair_id);im=ax.imshow(Z.median_effect.to_numpy()[:,None]*100,aspect='auto',cmap='coolwarm',vmin=-10,vmax=10);ax.set_xticks([0],[ring]);ax.set_yticks(range(len(core)),[r.niche+' '+r.celltype for r in core.itertuples()] if ring=='R1' else [],fontsize=7);ax.set_title(ring+' core-tier effects (pp)')
 fig.suptitle(scope+' | fixed core-tier overview; full 139 retained in index',fontsize=11);fig.colorbar(im,ax=list(axes),label='Median unit excess, pp');save(fig,root/'overview')
 log('PLOTS COMPLETE '+scope)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--scope',required=True,choices=SCOPES);args=ap.parse_args();plot_scope(args.scope)
