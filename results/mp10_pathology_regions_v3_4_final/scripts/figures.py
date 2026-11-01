"""Consistent panel coordinates and separate legend rows for the final atlas."""
from regions import *
from run_v3 import plt, PolyCollection, LineCollection, Patch
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'font.size':9})
LOCATION={'inside':('#bf3449','候选内侧'),'outside':('#4589b4','候选外侧'),'unlocalized':('#aaa4ab','位置未定'),'quality_excluded':('#383838','QC 排除')}
SUPPORT={'core':('#5b5697','作者可靠核心'),'parent':('#29978c','参考父区和组织邻接'),'enclosed':('#c87e2d','有完整包围证据的内部区室'),'conflict':('#cc6aa0','参考边缘仍有冲突'),'candidate':('#c7c0ad','沿用计算候选位置'),'weak':('#a38665','不稳健的小候选'),'unknown':('#999999','位置未定')}
def bands(w):
    return {
      'candidate_inside_beyond_w':('#a72740',f'候选内侧，距候选外缘 >{w} μm'),
      'inner_margin':('#eb966b',f'内侧 margin，≤{w} μm'),
      'outer_margin':('#ecc45e',f'外侧 margin，≤{w} μm'),
      'candidate_outside_beyond_w':('#4389b5',f'候选外侧，距候选外缘 >{w} μm'),
      'inside_outer_unobserved':('#c45865','内侧，未观察到外缘'),
      'outside_outer_unobserved':('#79a7c0','外侧，无可用外缘'),
      'outside_no_external_route':('#859db1','外侧区室，无可观测外侧路径'),
      'band_unresolved':('#e4e1dd','分带未定'),'unlocalized':('#999999','位置未定'),'quality_excluded':('#333333','QC 排除')}

def layout(rows,cols,title):
    fig=plt.figure(figsize=(4.5*cols,4.7*rows+0.7))
    gs=GridSpec(2*rows,cols,figure=fig,height_ratios=[5,1.45]*rows,left=.045,right=.975,bottom=.025,top=.91,hspace=.20,wspace=.22)
    axes=[];legends=[]
    for i in range(rows*cols):
        axes.append(fig.add_subplot(gs[(i//cols)*2,i%cols]))
        lg=fig.add_subplot(gs[(i//cols)*2+1,i%cols]);lg.axis('off');legends.append(lg)
    fig.suptitle(title,fontsize=12,y=.98)
    return fig,axes,legends

def legend(ax,handles):
    if handles:ax.legend(handles=handles,loc='center',ncol=2,fontsize=7,frameon=False,columnspacing=.9,handlelength=1.3)
def colored(ax,lg,polys,vals,palette):
    vals=np.asarray(vals)
    ax.add_collection(PolyCollection(polys,facecolors=[palette[v][0] for v in vals],edgecolors='none',antialiased=False))
    legend(lg,[Patch(color=c,label=l) for v,(c,l) in palette.items() if (vals==v).any()])
def setup(ax,title,limits,pitch,bar_um=500):
    ax.set_title(title,fontsize=10);ax.set_xlim(limits[0:2]);ax.set_ylim(limits[3],limits[2]);ax.axis('off');ax.set_aspect('equal')
    x=limits[0]+pitch*.5;y=limits[3]-pitch*.5
    ax.plot([x,x+pitch*bar_um/100],[y,y],color='#292929',lw=1.4)
    ax.text(x,y-pitch*.35,f'{bar_um} μm',fontsize=6,bbox={'facecolor':'white','edgecolor':'none','alpha':.7,'pad':.3})
def outlines(ax,e,all_types=False):
    lines=[('candidate_observed_outer','#252525','-')]
    if all_types:lines.extend([('sampling_or_tissue_cut_edge','#b37522','--'),('QC_or_missing_edge','#966e9e',':'),('internal_nonlesion_compartment_contact','#56847d',':')])
    for name,color,ls in lines:
        q=e[e.boundary_type==name]
        if len(q):
            v=q[['x1_px','y1_px','x2_px','y2_px']].to_numpy().reshape(-1,2,2)
            ax.add_collection(LineCollection(v,colors=color,linestyles=ls,linewidths=.65))
def save(fig,out,name):
    fig.savefig(out/(name+'.png'),dpi=300,facecolor='white');fig.savefig(out/(name+'.pdf'),facecolor='white');plt.close(fig)
def plot_section(sec,d,base,poly,image,e,co,out,ref,holes):
    pitch=co['median_neighbor_image_px']
    limits=[d.x.min()-2*pitch,d.x.max()+2*pitch,d.y.min()-2*pitch,d.y.max()+2*pitch]
    fig,a,lg=layout(2,2,sec+'  |  最终计算候选区域')
    for ax,title in zip(a,['最终位置','H&E 与候选外缘','距候选外缘的 200 μm 分带','位置支持与参考冲突']):setup(ax,title,limits,pitch)
    colored(a[0],lg[0],poly,d.lesion_location,LOCATION);outlines(a[0],e)
    a[1].imshow(image);outlines(a[1],e,True)
    legend(lg[1],[Line2D([0],[0],color=c,ls=s,label=t) for c,s,t in [('#252525','-','可观测候选外缘'),('#b37522','--','组织或捕获截断'),('#966e9e',':','QC 缺位'),('#56847d',':','内部区室接触')]])
    colored(a[2],lg[2],poly,d.distance_band_200um,bands(200));outlines(a[2],e)
    display=np.full(len(d),'candidate',object)
    display[d.position_support.eq('unstable_small_candidate')]='weak'
    display[d.reference_state.eq('possible_parent_reference_conflict')]='conflict'
    display[d.position_support.eq('author_parent_spatial_inference')]='parent'
    display[d.author_coarse_reference.isin(['author_coarse_inside','author_coarse_outside'])]='core'
    display[d.region_role.eq('lesion_internal_compartment')]='enclosed'
    display[~d.lesion_location.isin(['inside','outside'])]='unknown'
    colored(a[3],lg[3],poly,display,SUPPORT);save(fig,out,'final')
    fig,a,lg=layout(1,2,sec+'  |  同一候选轮廓的带宽对照')
    for ax,la,w in zip(a,lg,[100,300]):
        setup(ax,f'{w} μm',limits,pitch);colored(ax,la,poly,d[f'distance_band_{w}um'],bands(w));outlines(ax,e)
    save(fig,out,'widths')
    change=d.lesion_location.ne(base).to_numpy()
    if not change.any():return
    # Zoom around the specified objects, otherwise a changed connected patch.
    if sec.endswith('BC18'):
        groups=[ix for ix,pid in holes if len(ix)==225]
        ix=groups[0] if groups else np.flatnonzero(change)
    elif sec.endswith('BC13'):
        ix=np.flatnonzero(d.barcode.isin(['GATCCATAATCGCACG-1','GTACAGGACGTTACCT-1']))
    else:
        adj,_=graph(d);_,patches=components(adj,change);ix=max(patches,key=len)
    z=d.iloc[ix];span=max(12*pitch,z.x.max()-z.x.min()+6*pitch,z.y.max()-z.y.min()+6*pitch)
    cx=(z.x.min()+z.x.max())/2;cy=(z.y.min()+z.y.max())/2
    zoom=[cx-span*.6,cx+span*.6,cy-span*.6,cy+span*.6]
    fig,a,lg=layout(2,3,sec+'  |  参考父区域约束下的局部修订')
    for ax,title in zip([a[0],a[1]],['v3.3 位置','最终位置']):setup(ax,title,limits,pitch)
    colored(a[0],lg[0],poly,base,LOCATION);colored(a[1],lg[1],poly,d.lesion_location,LOCATION);outlines(a[1],e)
    a[2].imshow(ref[4]);a[2].axis('off');a[2].set_title('作者 S1 粗区域（作者坐标）',fontsize=10)
    legend(lg[2],[Patch(color=c,label=t) for c,t in [('#f59696','Tumor Cell'),('#5f88c9','Tumor Stroma'),('#5fadc9','Normal Stroma'),('#5fc9ad','Urothelium')]])
    for ax,title in zip(a[3:],['局部：v3.3','局部：最终位置','局部：最终 200 μm 分带']):setup(ax,title,zoom,pitch,200)
    colored(a[3],lg[3],poly,base,LOCATION);colored(a[4],lg[4],poly,d.lesion_location,LOCATION);colored(a[5],lg[5],poly,d.distance_band_200um,bands(200));outlines(a[4],e);outlines(a[5],e)
    save(fig,out,'revision')
