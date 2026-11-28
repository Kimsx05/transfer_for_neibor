from common3 import *
from scipy.signal import find_peaks
def main():
    p=read(O/'tables/section_predictions.tsv.gz');h=read(O/'tables/section_bins.tsv.gz')
    local=[]
    for (sec,b,dom,sg),a in p.groupby(['section_id','boundary','domain','segment']):
        a=a[a.x.abs()<=1000].sort_values('x')
        if len(a)<17 or a.method.iloc[0]!='penalized_cubic':continue
        for sign,what in [(1,'peak'),(-1,'valley')]:
            jj,info=find_peaks(sign*a.score,prominence=.005,distance=8)
            for j,pr,l,r in zip(jj,info['prominences'],info['left_bases'],info['right_bases']):
                if min(a.x.iloc[j]-a.x.iloc[l],a.x.iloc[r]-a.x.iloc[j])<100:continue
                local.append(dict(section_id=sec,boundary=b,domain=dom,segment=sg,kind=what,x=a.x.iloc[j],score=a.score.iloc[j],prominence=pr))
    ext=pd.DataFrame(local,columns=['section_id','boundary','domain','segment','kind','x','score','prominence']);save(ext,O/'tables/local_extrema.tsv')
    rows=[];cmp=[];status=read(O/'audit/section_boundary_status.tsv')
    for sec in M.index:
        for b in B:
            seg=p[(p.section_id==sec)&(p.boundary==b)&(p.domain=='E_inner')];a=seg[seg.x.abs()<=1000];ref=a[a.x==-50];found=ext[(ext.section_id==sec)&(ext.boundary==b)&(ext.domain=='E_inner')]
            connected=a[a.segment==ref.segment.iloc[0]] if len(ref) else a.iloc[:0]
            delta=connected.score.iloc[-1]-connected.score.iloc[0] if len(connected)>1 else np.nan
            hasmodel=a.method.isin(['penalized_cubic','constant']).any()
            kind='peak_and_valley' if found.kind.nunique()>1 else found.kind.iloc[0] if len(found) else 'no_curve' if not len(a) else 'bins_only' if not hasmodel else 'near_higher' if delta>=.005 else 'deep_higher' if delta<=-.005 else 'small_or_flat_change'
            rows.append(dict(section_id=sec,boundary=b,near_shape=kind,n_peaks=(found.kind=='peak').sum(),n_valleys=(found.kind=='valley').sum(),near_minus_farthest_contiguous=delta,valid_bins=h.query("section_id==@sec and boundary==@b and domain=='E_inner' and n>=5 and center>=-1000").shape[0],interface_status=status[(status.section_id==sec)&(status.boundary==b)].plot_status.iloc[0]))
        g=h[(h.section_id==sec)&h.n.ge(5)&h.center.abs().le(1000)]
        pair=g[g.domain=='E_inner'].pivot(index='center',columns='boundary',values='median').reindex(columns=B).dropna()
        out=g[(g.boundary==B[1])&g.center.gt(0)].pivot(index='center',columns='domain',values='median').reindex(columns=['EIN_EIS_outer','nonE_outer']).dropna()
        cmp.append(dict(section_id=sec,common_negative_bins=len(pair),rho=stats.spearmanr(pair.iloc[:,0],pair.iloc[:,1]).statistic if len(pair)>=3 and pair.nunique().min()>1 else np.nan,mean_abs_difference=pair.diff(axis=1).iloc[:,-1].abs().mean(),common_outer_bins=len(out),EIN_EIS_minus_nonE= (out.iloc[:,0]-out.iloc[:,1]).mean()))
    save(rows,O/'tables/section_shapes_near.tsv');save(cmp,O/'tables/shape_concordance.tsv')
if __name__=='__main__':main()

