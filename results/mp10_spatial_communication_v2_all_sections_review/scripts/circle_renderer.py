"""Common directed adjacency circle renderer. Pure display; never alters inference.
render_triplet(matrix, nodes, colors, epi_high, epi_low, out_prefix, title,
               subtitle='', computed=None, preview=False, metric='W')
matrix: source x target ndarray (NA retained); nodes includes both epi + targets.
Returns three metadata dictionaries. PNG 300 dpi and transparent vector PDF.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import textwrap

def render_triplet(matrix,nodes,colors,epi_high,epi_low,out_prefix,title,subtitle='',computed=None,preview=False,metric='W',epi_caption=None):
    mat=np.asarray(matrix,dtype=float).copy(); nodes=list(nodes); n=len(nodes)
    hi=nodes.index(epi_high);lo=nodes.index(epi_low);targets=[i for i in range(n) if i not in (hi,lo)]
    present=np.isfinite(mat).any(axis=0)|np.isfinite(mat).any(axis=1) if computed is None else np.asarray(computed,dtype=bool)
    # Fixed epi positions and target arc; shared by all three views including hidden-epi blank.
    xy=np.zeros((n,2));xy[hi]=[-.92,.63];xy[lo]=[-.92,-.63]
    theta=np.linspace(80,-80,len(targets)) if len(targets)>1 else np.array([0.])
    for k,t in enumerate(targets):xy[t]=[np.cos(np.deg2rad(theta[k])),np.sin(np.deg2rad(theta[k]))]
    keep=np.zeros((n,n),bool)
    for e in [hi,lo]:keep[e,targets]=True;keep[targets,e]=True
    masked=np.where(keep,mat,0.)
    vmax=float(np.nanmax(masked)) if np.isfinite(masked).any() else 0.;vmax=max(0.,vmax)
    results=[];prefix=Path(out_prefix);prefix.parent.mkdir(parents=True,exist_ok=True)
    width=12.;height=max(9.,.7*len(targets)+2.5)
    for view,states in [('High',[hi]),('Low',[lo]),('Both',[hi,lo])]:
        fig,ax=plt.subplots(figsize=(width,height));fig.subplots_adjust(left=.04,right=.97,bottom=.14,top=.86)
        active=targets+states;positive=[]
        for e in states:
            for t in targets:
                for a,b in [(e,t),(t,e)]:
                    value=mat[a,b]
                    if np.isfinite(value) and value>0:
                        # Opposite directions use the same signed curvature relative to their
                        # own orientation, hence lie on opposite sides of the straight chord.
                        lw=6*value/vmax
                        edge=FancyArrowPatch(xy[a],xy[b],connectionstyle='arc3,rad=0.13',arrowstyle='-|>',mutation_scale=12,linewidth=lw,color=colors[a],alpha=.70,shrinkA=13,shrinkB=13,zorder=1)
                        ax.add_patch(edge);positive.append((a,b,value))
        for i in active:
            ax.scatter(*xy[i],s=310,facecolor=colors[i] if present[i] else '#eeeeee',edgecolor='black' if present[i] else '#999999',linewidth=.9,zorder=3)
            label=nodes[i]+(' [NE]' if not present[i] else '')
            if i in (hi,lo):label='Epi_H' if i==hi else 'Epi_L'
            label='\n'.join(textwrap.wrap(label,width=31))
            if i in (hi,lo):tx,ty=xy[i]+[-.09,0];ha='right'
            else:tx,ty=xy[i]+[.09,0];ha='left'
            ax.text(tx,ty,label,ha=ha,va='center',fontsize=10.5,fontweight='bold' if nodes[i]=='Treg' else 'normal',color='black',zorder=4)
        ax.set_xlim(-1.70,2.13);ax.set_ylim(-1.19,1.19);ax.axis('off')
        fig.suptitle(f'{title} | {view}\nTotal interaction strength',fontsize=15,y=.96,color='black',ha='center')
        # Exact linear mapping (no visibility floor), with explicit values; no edge selection.
        legend_ax=fig.add_axes([.33,.073,.48,.023]);legend_ax.set_xlim(0,1);legend_ax.set_ylim(0,1);legend_ax.axis('off')
        for x,f in zip([.04,.37,.70],[.25,.5,1.]):
            legend_ax.plot([x,x+.17],[.72,.72],lw=6*f,color='#555555');legend_ax.text(x+.085,-.12,f'{vmax*f:.4g}',ha='center',va='center',fontsize=9)
        fig.text(.5,.108,f'{metric}: shared High/Low/Both linear width scale; max = {vmax:.6g}',ha='center',fontsize=10,color='black')
        fig.text(.5,.043,epi_caption or 'Epi_H / Epi_L = MP10-high / low epithelial states | arrows: source → target | equal node sizes',ha='center',fontsize=9,color='black')
        fig.text(.5,.019,subtitle or 'All computed LR, internal P < 0.05; NE is unavailable, not zero. Scales differ across panels.',ha='center',fontsize=8.5,color='black',wrap=True,va='center' if '\n' in subtitle else 'baseline')
        stem=Path(str(prefix)+'_'+view);fig.savefig(str(stem)+'.pdf',transparent=True);fig.savefig(str(stem)+'.png',transparent=True,dpi=300)
        if preview and view=='Both':fig.savefig(str(stem)+'_white.png',facecolor='white',dpi=130)
        plt.close(fig)
        results.append(dict(view=view,pdf=str(stem)+'.pdf',png=str(stem)+'.png',edge_max=vmax,n_positive_edges=len(positive),n_nodes=len(active),n_computed_nodes=sum(present[active]),metric=metric,hidden_positive_edges=0,has_self_loops=False))
    return results
