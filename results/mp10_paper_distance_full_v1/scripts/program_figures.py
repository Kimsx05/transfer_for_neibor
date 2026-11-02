from figures import *
def maps(s):
    d=D[s]
    for page,cols in enumerate([focus[:6],focus[6:12],["Treg_fraction","Treg_abundance",focus[-2],focus[-1],focus[11],focus[6]]]):
        title=["Growth, stress and inflammation","Immune, lipid and matrix programs","Treg context across the full tissue"][page]
        f,a,l=layout(2,3,s+"  |  "+title,"Whole-spot RNA programs (UCell)   •   display scales fixed per feature across all sections")
        for j,c in enumerate(cols):
            unit="UCell score";title=short(c)
            if c=="Treg_abundance":unit="Posterior mean estimate";title="Treg estimated abundance"
            if c=="Treg_fraction":unit="Fraction of all 30";title="Treg composition"
            field(a[j],l[j],d,c,title,unit)
        export(f,a,l,s,f"0{5+page}_program_maps")
def trends(s):
    d=D[s]; modules=[focus[i] for i in [4,2,0,1,9,10]]
    f,a,l=layout(2,3,s+"  |  Programs by target distance","Epi queries   •   medians at observed distance bins   •   lines stop where a bin has fewer than 5 spots",h=3,w=4.7)
    for j,c in enumerate(modules):
        ax=a[j]
        for i,t in enumerate(targets):
            x=d[t+"_P80_all_distance_um"].round(6)
            b=pd.cut(x,[-1,0,200,500,1000,np.inf],labels=False)
            g=d[c].groupby(b).agg(["median","count"]).reindex(range(5))
            ax.plot(range(5),g["median"].where(g["count"]>=5),color=colors_target[i],marker="o",ms=4,label=t)
        ax.set_title(short(c));ax.set_xticks(range(5),["0","≤200","≤500","≤1000",">1000"],rotation=15)
        ax.set_xlabel("Distance-bin upper bound (μm)");ax.set_ylabel("Median spot UCell")
        l[j].legend(handles=ax.get_lines(),loc="center",ncol=3,frameon=False,handlelength=1)
    export(f,a,l,s,"08_program_trends")
    z=pd.read_csv(tab/"program_distance.tsv",sep="\t")
    z=z[z.section_id==s].pivot(index="program",columns="target",values="pearson_r")
    f,a,l=layout(1,2,s+"  |  Program–distance associations","All 55 UCell programs   •   one section   •   including zero distances   •   display cap ±0.6",w=8,h=12,left=.27,wspace=1.85)
    for j,mat in enumerate([z.iloc[:28],z.iloc[28:]]):
        mat=mat.copy();mat.index=[short(x) for x in mat.index]
        heat(a[j],l[j],mat,"Whole-spot programs in Epi queries",-.6,.6,"RdBu_r","Pearson r")
    export(f,a,l,s,"09_program_effects")
