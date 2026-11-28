from common3 import *
from scipy.interpolate import BSpline
def smooth(x,y):
    # Uniform cubic B-splines, second-difference penalty one, fixed model size.
    n=len(x);k=min(5,n-1);t=np.r_[[x[0]]*4,np.linspace(x[0],x[-1],k-2)[1:-1],[x[-1]]*4]
    design=BSpline.design_matrix(x,t,3).toarray();pen=np.diff(np.eye(k),n=2,axis=0)
    coef=np.linalg.solve(design.T@design+pen.T@pen,design.T@y)
    xx=np.linspace(x[0],x[-1],(n-1)*8+1);return xx,BSpline(t,coef,3)(xx)
def prepare():
    bins=[];pred=[];coverage=[];distributions=[];contrib=[];cache=[]
    cols=['MP10','new_epi_domain','expression_QC','region_main','all_E',B[1],'P95_new','E_f','T_f','I_f','UMI','block','legacy_query']+[k+b for b in B for k in ['dist_','strict_']]+[b+'_component' for b in B]
    for sec in M.index:
        d=spot(sec,cols);v=d.expression_QC&d.MP10.notna();q=v&d.all_E;assert q.equals(d.new_epi_domain)
        mu=d.loc[q,'MP10'].mean();cache.append(dict(section_id=sec,spots_path=str(V/'spots'/(sec+'.parquet')),boundary_path=str(V/'boundary'/(sec+'.parquet')),n=len(d),expression_QC=int(d.expression_QC.sum()),finite=int(d.MP10.notna().sum()),D_MP10=int(q.sum()),old_hidden=int((q&~d.legacy_query).sum()),zero=int((v&d.MP10.eq(0)).sum()),center=mu,raw_max=d.MP10[v].max(),center_absmax=(d.MP10[v]-mu).abs().max()))
        for b in B:
            dist=d['dist_'+b];ok=v&dist.notna();inside=d[b];assert not ((dist<0)&~inside|(dist>0)&inside).any()
            assert not ((dist>0)&q&(~d.region_main.isin(['EIN','EIS']) if b==B[1] else True)).any()
            labels=np.where(inside,'E_inner',np.where(d.all_E,'EIN_EIS_outer','nonE_outer'))
            k=dict(section_id=sec,boundary=b)
            coverage.append(dict(k,mask_n=inside.sum(),n_total=len(d),qc_fail=(~d.expression_QC).sum(),no_score=(~np.isfinite(d.MP10)).sum(),unknown=(d.region_main=='Unknown').sum(),distance_NA=dist.isna().sum(),display_n=ok.sum(),plot_status='ok' if ok.any() else 'empty_mask' if not inside.any() else 'no_usable_observed_interface',zero_n=(ok&d.MP10.eq(0)).sum()))
            ix=np.maximum(0,np.ceil(np.abs(dist[ok].to_numpy())/100-1e-10)-1)
            a=d.loc[ok].assign(domain=labels[ok],center=np.where(inside[ok],-1,1)*(ix*100+50),logUMI=np.log1p(d.UMI[ok]))
            for (dom,x),g in a.groupby(['domain','center']):
                j=dict(k,domain=dom,center=x);z=g.MP10
                bins.append(dict(j,n=len(g),median=z.median(),q25=z.quantile(.25),q75=z.quantile(.75),mean=z.mean(),p95_fraction=g.P95_new.mean() if dom!='nonE_outer' else np.nan,centered=z.median()-mu,baseline=mu,E_fraction=g.E_f.mean(),T_fraction=g.T_f.mean(),I_fraction=g.I_f.mean(),logUMI=g.logUMI.mean(),n_blocks=g.block.nunique(),n_E_components=g.loc[g[b],b+'_component'].nunique()))
                for (block,comp),nn in g.groupby(['block',b+'_component']).size().items():contrib.append(dict(j,block=block,component=comp,n=nn))
            for dom,g in a.groupby('domain'):
                distributions.append(dict(k,domain=dom,n=len(g),zero=(g.MP10==0).sum(),mean=g.MP10.mean(),median=g.MP10.median(),near_n=(g.center.abs()<=1000).sum()))
    h=pd.DataFrame(bins);shapes=[]
    for (sec,b,dom),a in h[h.n>=5].groupby(['section_id','boundary','domain']):
        a=a.sort_values('center');mu=a.baseline.iloc[0]
        for seg,g in a.groupby(a.center.diff().ne(100).cumsum()):
            x=g.center.to_numpy();y=g['median'].to_numpy();method='bin_summary';xx=x;yy=y
            if len(y)>=5:
                if np.ptp(y)==0:method='constant'
                else:xx,yy=smooth(x,y);method='penalized_cubic'
            dy=np.diff(yy);directions=np.sign(dy[np.abs(dy)>1e-10]);turn=np.diff(directions);j=dict(section_id=sec,boundary=b,domain=dom,segment=int(seg))
            shapes.append(dict(j,n_bins=len(x),xmin=x.min(),xmax=x.max(),method=method,range=float(np.ptp(yy)),peaks=sum(turn<0),valleys=sum(turn>0),max_at=xx[np.argmax(yy)],min_at=xx[np.argmin(yy)],edge_minus_far=(yy[-1]-yy[0])*(1 if dom=='E_inner' else -1)))
            pred.extend(dict(j,x=u,score=v,centered=v-mu,method=method) for u,v in zip(xx,yy))
    save(cache,O/'cache/point_index.tsv');save(h,O/'tables/section_bins.tsv.gz');save(pred,O/'tables/section_predictions.tsv.gz');save(shapes,O/'tables/section_shapes.tsv.gz');save(contrib,O/'tables/bin_geometry_support.tsv.gz');save(distributions,O/'tables/display_domains.tsv');save(coverage,O/'audit/section_boundary_status.tsv')
def plot():
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    theme();s=read(O/'cache/point_index.tsv').set_index('section_id');h=read(O/'tables/section_bins.tsv.gz');p=read(O/'tables/section_predictions.tsv.gz');status=read(O/'audit/section_boundary_status.tsv')
    path=O/'figures/sections';path.mkdir(exist_ok=True,parents=True);rows=[];pal={'nonE_outer':'#888989','EIN_EIS_outer':'#aa6186','E_inner':'#267689'}
    for sec in M.index:
        d=spot(sec,['MP10','expression_QC','all_E',B[1]]+['dist_'+b for b in B]);mu=s.at[sec,'center']
        for b in B:
            st=status[(status.section_id==sec)&(status.boundary==b)].iloc[0];z=d[d.expression_QC&d.MP10.notna()&d['dist_'+b].notna()].copy();di=z['dist_'+b];z['domain']=np.where(z[b],'E_inner',np.where(z.all_E,'EIN_EIS_outer','nonE_outer'));a=h[(h.section_id==sec)&(h.boundary==b)&(h.n>=5)];pp=p[(p.section_id==sec)&(p.boundary==b)]
            for mode,view in itertools.product(['raw','centered'],['near','full']):
                key=dict(section_id=sec,boundary=b,scale=mode,view=view,status=st.plot_status,n=len(z))
                if not len(z):rows.append(dict(key,file=None));continue
                target=path/f'{sec}__{b}__{mode}_{view}.pdf';rows.append(dict(key,file=str(target.relative_to(O))))
                off=0 if mode=='raw' else mu
                fig,ax=plt.subplots(figsize=(6.7,4.45));fig.subplots_adjust(left=.13,right=.945,bottom=.17,top=.89)
                for dom,col in pal.items():
                    g=z[z.domain==dom];t=a[a.domain==dom];pr=pp[pp.domain==dom]
                    if not len(g):continue
                    ax.scatter(g['dist_'+b],g.MP10-off,s=4,alpha=.23,lw=0,c=col,rasterized=True)
                    ax.scatter(t.center,t['median']-off,s=7,c=col,zorder=5)
                    for _,ss in pr.groupby('segment'):
                        hh=t[(t.center>=ss.x.min())&(t.center<=ss.x.max())].sort_values('center')
                        ax.fill_between(hh.center,hh.q25-off,hh.q75-off,color=col,alpha=.15);ax.plot(ss.x,ss.score-off,color=col,lw=1.45)
                    ax.plot([],[],color=col,label={'E_inner':'E interior','nonE_outer':'non-E background','EIN_EIS_outer':'EIN/EIS exterior'}[dom])
                ax.axvline(0,color='#303030',lw=.8,ls=':');ax.set_xlim((1000,-1000) if view=='near' else (max(0,di.max())+40,min(0,di.min())-40))
                lim=s.raw_max.max()*1.02 if mode=='raw' else s.center_absmax.max()*1.02;ax.set_ylim(0 if mode=='raw' else -lim,lim)
                ax.set(xlabel='Signed distance (μm)',ylabel='MP10 score' if mode=='raw' else 'Centered MP10 score');ax.set_title(sec+' · '+('all E' if b==B[0] else 'ENN + ENS'),fontsize=10);handles,_=ax.get_legend_handles_labels();handles.append(Patch(facecolor='#c5c5c5',alpha=.55,label='Band: IQR'));ax.legend(handles=handles,loc='upper right',fontsize=7,frameon=False);fig.savefig(target,dpi=125);plt.close(fig)
    save(rows,O/'audit/section_figure_index.tsv')
if __name__=='__main__':prepare() if 'prepare' in sys.argv else plot()

