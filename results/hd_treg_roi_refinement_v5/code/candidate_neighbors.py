"""One observed physical neighborhood per RNA-supported candidate, all ten, ordered by ID."""
from pictures import View
from common5 import R, pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors

s = 'D1_Johnson'
view = View(s)
d = pd.read_parquet(R/s/'cells_scored.parquet')
ne = pd.read_parquet(R/s/'neighborhoods.parquet')
cases = d[d.group.eq('RNA_candidate')].sort_values('id')
epi = d[d.epi_eligible]
norm = colors.Normalize(0, epi.MP10_cell_score.max())
with PdfPages(R/'figures/D1_candidate_neighborhoods.pdf') as pdf:
    for start in range(0, len(cases), 5):
        fig, axs = plt.subplots(2,5,figsize=(16,7),layout='constrained',squeeze=False)
        for j,z in enumerate(cases.iloc[start:start+5].itertuples()):
            x,y = z.x,z.y
            ax = axs[0,j]
            ax.imshow(view.sample(x,y,60,240,'edge'),extent=[-30,30,30,-30])
            ax.set_title(f'{z.id}  {z.pred}',fontsize=12)
            ax.plot([13,23],[24,24],color='white',linewidth=2)
            ax.set_axis_off()
            ax=axs[1,j];extent=220
            ax.imshow(view.sample(x,y,2*extent,440,'he'),extent=[-extent,extent,extent,-extent],alpha=.33)
            g=epi[epi.component.eq(z.component)&epi.x.sub(x).abs().le(extent)&epi.y.sub(y).abs().le(extent)]
            ok=g[g.score_ok];u=g[~g.score_ok]
            ax.scatter(u.x-x,u.y-y,marker='x',s=8,color='grey',linewidths=.5)
            pts=ax.scatter(ok.x-x,ok.y-y,c=ok.MP10_cell_score,s=9,cmap='viridis',norm=norm,linewidths=0)
            h=ok[ok.H05]
            ax.scatter(h.x-x,h.y-y,s=22,edgecolors='#c55020',facecolors='none',linewidths=1)
            ax.scatter([0],[0],c='#0077bb',marker='+',s=60,linewidths=1.3)
            for rad in (50,200):ax.add_patch(Circle((0,0),rad,fill=False,color='#333333',linewidth=1))
            n=ne[ne.id.eq(z.id)].sort_values('radius')
            ax.set_xlabel('50 / 200 µm scored epithelium\n'+' / '.join(map(str,n.scored))+', coverage '+' / '.join(f'{w:.0%}' for w in n.coverage),fontsize=9)
            ax.set_xticks([]);ax.set_yticks([])
        fig.colorbar(pts,ax=list(axs[1,:]),location='bottom',shrink=.5,pad=.07,label='Epithelial cell MP10   |   orange ring = H05   |   + = RNA candidate')
        pdf.savefig(fig,dpi=220)
        fig.savefig(R/'figures'/f'D1_candidate_neighborhoods_{start//5+1}.png',dpi=160)
        plt.close(fig)
