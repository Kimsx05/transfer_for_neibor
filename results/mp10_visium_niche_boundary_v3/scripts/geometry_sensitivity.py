from common3 import *
from down import locked,unit_bins,runmodel
def build():
    out=[]
    for s in M.index:
        d=spot(s);q=d.new_epi_domain;mu=d.MP10[q].mean();common=q&d[B[1]]&d.dist_all_E.notna()&d['dist_'+B[1]].notna()
        for b in B:
            distance=d['dist_'+b];near=d[b]&q&distance.ge(-200)&distance.le(0);deep=d[b]&q&distance.lt(-200);comp=d[b+'_component']
            ns=comp[near].value_counts();ds=comp[deep].value_counts();keep=ns[ns>=3].index.intersection(ds[ds>=3].index)
            for variant,dom,di in [('common_points',common,distance),('strict',q&d[b],d['strict_'+b]),('two_band_components',q&d[b]&comp.isin(keep),distance)]:
                hh=d[dom&di.notna()].copy();hh['center']=-(np.maximum(0,np.ceil(np.abs(di.loc[hh.index])/100-1e-10)-1)*100+50);hh['logUMI']=np.log1p(hh.UMI)
                for x,v in hh.groupby('center'):out.append(dict(section_id=s,boundary=b,domain='E_inner',variant=variant,center=x,n=len(v),median=v.MP10.median(),centered=v.MP10.median()-mu,E_fraction=v.E_f.mean(),logUMI=v.logUMI.mean()))
    save(out,O/'cache/geometry_bins.tsv.gz')
def analyze():
    lab=locked();t=read(O/'cache/geometry_bins.tsv.gz');sup=[];test=[];diff=[];groups=[]
    for var,a in t.groupby('variant'):
        aa,u=unit_bins(a,lab)
        for b in B:
            z=u[u.boundary==b];orig=aa[aa.boundary==b]
            runmodel(z,orig,b,'E_inner',G,var,sup,test,diff)
            for keys,g in z.groupby(['group','center']):
                x=g.groupby('dataset_id')[['raw','centered']].mean().mean()
                groups.append(dict(variant=var,boundary=b,group=keys[0],center=keys[1],raw=x.raw,centered=x.centered,n_units=g.unit.nunique(),n_sources=g.dataset_id.nunique()))
    r=pd.DataFrame(test);r['q']=np.nan
    for _,i in r.groupby(['variant','term']).groups.items():r.loc[i,'q']=bh(r.loc[i,'p'])
    save(r,OUT/'tables/geometry_curve_tests.tsv');save(sup,OUT/'tables/geometry_curve_support.tsv.gz');save(diff,OUT/'tables/geometry_curve_differences.tsv.gz');save(groups,OUT/'tables/geometry_group_trends.tsv.gz')
if __name__=='__main__':build() if 'prepare' in sys.argv else analyze()

