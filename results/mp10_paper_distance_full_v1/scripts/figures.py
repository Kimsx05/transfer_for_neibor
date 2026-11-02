"""Generate white vector PDFs with separate bodies/keys, and 300-dpi PNGs."""
from distance import *
import gzip,json,math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import LineCollection
from matplotlib.transforms import Bbox
from PIL import Image
from matplotlib import font_manager
import subprocess
for style in ["Regular","Bold"]:
    font_manager.fontManager.addfont(subprocess.check_output(["fc-match","-f","%{file}","Arial:style="+style],text=True))
plt.rcParams.update({"font.family":"Arial","axes.linewidth":1.1,"font.size":12,"axes.titlesize":12,"axes.labelsize":12,"legend.fontsize":11,"pdf.fonttype":42,"lines.linewidth":1.1,"figure.facecolor":"white","axes.facecolor":"white","axes.spines.right":False,"axes.spines.top":False})
manifest=pd.read_csv(out/"tables/input_manifest.tsv",sep="\t")
secs=manifest.section_id.tolist()
D={}
scales={}
section_ds=manifest.set_index("section_id").dataset_id.to_dict()
A=pd.read_csv(tab/"associations.tsv",sep="\t")
R=pd.read_csv(tab/"region_composition.tsv",sep="\t")
N=pd.read_csv(tab/"region_niche.tsv",sep="\t")
P=pd.read_csv(root/"neighborhood_niche_visualization_v1/data_tables/niche_palette.tsv",sep="\t").set_index("accepted_label").color.to_dict()
regions={"inside":"#bf3449","outside":"#4589b4","unlocalized":"#aaa4ab","quality_excluded":"#383838"}
colors_target=["#367597","#A85417","#BD4260"]
focus=["HALLMARK_"+i for i in ["E2F_TARGETS","G2M_CHECKPOINT","HYPOXIA","TNFA_SIGNALING_VIA_NFKB","EPITHELIAL_MESENCHYMAL_TRANSITION","INTERFERON_GAMMA_RESPONSE","IL2_STAT5_SIGNALING","TGF_BETA_SIGNALING","CHOLESTEROL_HOMEOSTASIS"]]+["REACTOME_"+i for i in ["CHOLESTEROL_BIOSYNTHESIS","EXTRACELLULAR_MATRIX_ORGANIZATION","CO_INHIBITION_BY_CTLA4"]]+["GSE15659_RESTING_VS_ACTIVATED_TREG_DN","GSE15659_RESTING_VS_ACTIVATED_TREG_UP"]
def short(s):
    return s.replace("HALLMARK_","").replace("REACTOME_","").replace("GSE15659_RESTING_VS_ACTIVATED_TREG_","Treg resting vs activated ").replace("EPITHELIAL_MESENCHYMAL_TRANSITION","EMT").replace("_"," ").title()
def layout(rows,cols,title,subtitle,w=4.3,h=4,left=.065,wspace=.28):
    fig=plt.figure(figsize=(cols*w,rows*(h+.85)+1.1))
    gs=fig.add_gridspec(2*rows,cols,height_ratios=[h,.72]*rows,left=left,right=.975,top=1-1.35/fig.get_figheight(),bottom=.035,hspace=.55,wspace=wspace)
    fig.text(.065,1-.15/fig.get_figheight(),title,fontsize=16,fontweight="bold",va="top")
    fig.text(.065,1-.54/fig.get_figheight(),subtitle,fontsize=12,va="top")
    ax=[];lg=[]
    for i in range(rows*cols):
        ax.append(fig.add_subplot(gs[2*(i//cols),i%cols]))
        lg.append(fig.add_subplot(gs[2*(i//cols)+1,i%cols]));lg[-1].axis("off")
    return fig,ax,lg
def export(fig,axes,keys,s,name):
    dest=out/"figures"/section_ds[s]/s if s in section_ds else out/"figures"/s;dest.mkdir(parents=True,exist_ok=True)
    path=dest/name
    fig.savefig(str(path)+".pdf",facecolor="white")
    fig.savefig(str(path)+".png",dpi=300,facecolor="white")
    for a in keys:a.set_visible(False)
    fig.savefig(str(path)+"_panel_only.pdf",facecolor="white",bbox_inches="tight")
    for a in keys:a.set_visible(True)
    for a in axes:a.set_visible(False)
    for t in fig.texts:t.set_visible(False)
    fig.savefig(str(path)+"_legend_only.pdf",facecolor="white",bbox_inches="tight")
    plt.close(fig)
def keybar(ax,cmap,lo,hi,label):
    # Draw a vector strip in a dedicated legend axis.
    pal=plt.get_cmap(cmap)
    for j in range(64):ax.add_patch(plt.Rectangle((.08+.84*j/64,.59),.84/64+.0001,.18,transform=ax.transAxes,color=pal((j+.5)/64),linewidth=0))
    for x,v in [(.08,lo),(.92,hi)]:ax.text(x,.34,f"{v:.3g}",transform=ax.transAxes,ha="center",va="center",fontsize=11)
    ax.text(.5,.02,label,transform=ax.transAxes,ha="center",va="center",fontsize=11)
def category(ax,pal):
    if not pal:
        ax.text(.5,.5,"No covered labels",ha="center",transform=ax.transAxes);return
    ax.legend(handles=[Patch(facecolor=c,label=k) for k,c in pal.items()],loc="center",ncol=min(5,len(pal)),frameon=False,columnspacing=.7,handlelength=1)
def setup(ax,d,title):
    ax.set_title(title)
    ax.set_xlim(d.x.min()-5,d.x.max()+5);ax.set_ylim(d.y.max()+5,d.y.min()-5);ax.set_aspect("equal");ax.axis("off")
def bounds(col):
    return (0,scales[col])
def field(ax,lg,d,col,title,unit,limits=None,cmap="viridis",mask=None):
    setup(ax,d,title);ax.scatter(d.x,d.y,s=7,color="#dedede",linewidths=0)
    u=d[col].notna()
    if mask is not None:u &= mask
    lo,hi=limits or bounds(col)
    ax.scatter(d.x[u],d.y[u],s=7,c=d.loc[u,col],cmap=cmap,vmin=lo,vmax=hi,linewidths=0)
    keybar(lg,cmap,lo,hi,unit)
def catmap(ax,lg,d,col,pal,title):
    setup(ax,d,title);ax.scatter(d.x,d.y,s=7,c=[pal.get(v,"#dedede") for v in d[col].astype(object)],linewidths=0)
    category(lg,{str(k):v for k,v in pal.items() if (d[col]==k).any()})
def core(sec):
    d=D[sec];r=manifest.set_index("section_id").loc[sec]
    f,a,l=layout(2,3,sec+"  |  Tissue context","Aligned image pixels for display   •   frozen v3.4 regions   •   grey = outside displayed coverage")
    setup(a[0],d,"H&E and observed lesion edge")
    a[0].imshow(Image.open(r.image_path))
    path=root/"results/mp10_pathology_regions_v3_4_final/by_dataset"/r.dataset_id/sec/"geometry.json.gz"
    with gzip.open(path,"rt") as src:geom=json.load(src)
    e=pd.DataFrame(geom["roi_edges"]["data"],columns=geom["roi_edges"]["columns"])
    e=e[e.boundary_type=="candidate_observed_outer"]
    a[0].add_collection(LineCollection(e[["x1_px","y1_px","x2_px","y2_px"]].to_numpy().reshape(-1,2,2),colors="#222222",linewidths=1))
    l[0].text(.5,.5,"Black line: observed candidate outer edge" if len(e) else "No observed candidate outer edge",ha="center",transform=l[0].transAxes)
    catmap(a[1],l[1],d,"lesion_location",regions,"Final candidate location")
    catmap(a[2],l[2],d,"accepted_label",P,"Accepted R12 (original 0–8)")
    field(a[3],l[3],d,"Epithelial_fraction","Epithelial support","Fraction of all 30",(0,1))
    field(a[4],l[4],d,raw,"MP10 in Epi queries","Mixed-spot Top20 AUCell",mask=d["query"])
    field(a[5],l[5],d,"epi_MP10_model_universe","Epi-component MP10 (valid only)","Matched-model AUCell",mask=d.epi_score_valid.eq(True))
    export(f,a,l,sec,"01_context")
    f,a,l=layout(3,3,sec+"  |  Candidate neighborhoods","Columns: full-tissue composition, P80 enriched locations, nearest target distance in Epi queries")
    for i,t in enumerate(targets):
        j=i*3;tag=t+"_P80_all"
        field(a[j],l[j],d,t+"_fraction",t+" composition","Fraction of all 30")
        catmap(a[j+1],l[j+1],d,tag+"_target",{False:"#dcdddf",True:colors_target[i]},t+" P80 targets")
        field(a[j+2],l[j+2],d,tag+"_distance_um",t+" distance","Nearest spot center (μm)")
    export(f,a,l,sec,"02_targets")
    f,a,l=layout(3,3,sec+"  |  MP10 versus nearest target","Mixed-spot Top20 AUCell in Epi queries   •   section P80, all tissue targets   •   same tissue component",w=4.8,h=3.4)
    for i,t in enumerate(targets):
        x=d[t+"_P80_all_distance_um"];v=d[raw]
        stats=A[(A.section_id==sec)&(A.target==t)&(A["quantile"]==80)&(A.scope=="all")]
        for positive in [False,True]:
            j=3*i+int(positive);u=x.gt(0) if positive else x.notna()
            ax=a[j];ax.scatter(x[u],v[u],s=7,c=np.where(x[u]==0,"#d55e00",colors_target[i]),linewidths=0,alpha=.48)
            if x[u].nunique()>1:
                xx=np.array([x[u].min(),x[u].max()]);ax.plot(xx,np.polyval(np.polyfit(x[u],v[u],1),xx),c="#222222",lw=1.6)
            ax.set_title(t+("  |  positive distances" if positive else "  |  including overlap"));ax.set_xlabel("Nearest target (μm)");ax.set_ylabel("MP10 AUCell")
            z=stats[stats.positive_only==positive].iloc[0]
            l[j].text(.5,.5,f"r {z.pearson_r:+.2f}    ρ {z.spearman_rho:+.2f}\nn {z.n_used}/{z.n_query} queries   •   {z.n_target} targets",transform=l[j].transAxes,ha="center",va="center")
        ax=a[3*i+2];dist=x.round(6);bb=pd.cut(dist,[-1,0,200,500,1000,np.inf],labels=["0","≤200","≤500","≤1000",">1000"])
        groups=[v[bb.eq(b)] for b in bb.cat.categories]
        ax.boxplot(groups,showfliers=False,patch_artist=True,boxprops={"facecolor":colors_target[i]},medianprops={"color":"black"})
        ax.set_xticklabels([f"{b}\n{len(g)}" for b,g in zip(bb.cat.categories,groups)],fontsize=10);ax.set_title(t+"  |  distance bins");ax.set_xlabel("Upper bound (μm) / spot count");ax.set_ylabel("MP10 AUCell")
        l[3*i+2].text(.5,.5,"Orange = zero distance\nLine = within-section linear fit",transform=l[3*i+2].transAxes,ha="center",va="center")
    export(f,a,l,sec,"03_distance")
    region(sec)
def heat(ax,lg,mat,title,lo,hi,cmap,unit):
    ax.pcolormesh(np.ma.masked_invalid(mat.to_numpy(dtype=float)),cmap=cmap,vmin=lo,vmax=hi,shading="flat")
    ax.set_facecolor("#ddd");ax.set_xlim(0,len(mat.columns));ax.set_ylim(len(mat),0)
    ax.set_xticks(np.arange(len(mat.columns))+.5,mat.columns,rotation=48,ha="right",fontsize=10)
    ax.set_yticks(np.arange(len(mat))+.5,mat.index,fontsize=10);ax.set_title(title);keybar(lg,cmap,lo,hi,unit)
def region(sec):
    f,a,l=layout(2,2,sec+"  |  Region composition","v3.4 200 μm bands   •   each mean gives equal weight to observed spots   •   grey = unavailable",w=8.3,h=4.1,left=.16,wspace=.64)
    z=R[(R.section_id==sec)&(R.axis=="distance_band_200um")].set_index("region")
    z=z.reindex([k for k in region_order if k in z.index])
    names={"candidate_inside_beyond_w":"Inside >200","candidate_outside_beyond_w":"Outside >200","band_unresolved":"Unresolved","outside_no_external_route":"Outside, no route","inside_outer_unobserved":"Inside, no edge","outside_outer_unobserved":"Outside, no edge"}
    def rename(df):return df.rename(index=names).rename(columns=lambda x:x.replace("_fraction","").replace("_"," "))
    majors=[c for c in z if c.endswith("_fraction") and c.split("_fraction")[0] not in targets]
    heat(a[0],l[0],rename(z[majors]),"Reference lineages",0,1,"viridis","Mean fraction of all 30")
    heat(a[1],l[1],rename(z[[t+"_fraction" for t in targets]+["Treg_over_T_NK"]]),"Candidates and Treg within T/NK",0,.35,"magma","Mean spot fraction (display cap 0.35)")
    z=N[(N.section_id==sec)&(N.version=="SAMPLE_R12_LEIDEN")&(N.axis=="distance_band_200um")]
    for j,col,title,limit in [(2,"fraction","R12 among labelled Epi anchors",1),(3,"MP10_median","MP10 by band and R12 (n ≥ 5)",bounds(raw)[1])]:
        mat=z.pivot(index="region",columns="niche",values=col)
        mat=mat.reindex([k for k in region_order if k in mat.index])
        n=z.groupby("region").first()
        mat.index=[names.get(x,x)+f" ({n.loc[x,'n_labelled']}/{n.loc[x,'n_eligible']})" for x in mat.index]
        heat(a[j],l[j],mat,title,0,limit,"viridis",("Composition" if j==2 else "Mixed-spot AUCell"))
    export(f,a,l,sec,"04_regions")
