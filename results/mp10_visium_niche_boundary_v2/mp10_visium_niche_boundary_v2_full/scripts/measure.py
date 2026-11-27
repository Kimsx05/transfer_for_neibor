from lib import *
assert (O/'classification/classification_lock.json').is_file()
effects=[];counts=[];bins=[];geos=[];overall=[];areas=[];comps=[];paired=[];outer=[];clusters=[]
def contrast(d,domain,b,dist,w,variant,inside):
    p=np.asarray(domain)&np.asarray(inside)&np.isfinite(dist)
    n=p&(np.abs(dist)<=w);deep=p&(np.abs(dist)>w);n1=n.sum();n2=deep.sum()
    f=np.asarray(d.MP10);value=avg(f[n])-avg(f[deep])
    pos=d.P95_new.to_numpy(dtype=float,na_value=np.nan);prob=100*(avg(pos[n&np.isfinite(pos)])-avg(pos[deep&np.isfinite(pos)]))
    n_p95_B1=sum(n&np.isfinite(pos));n_p95_B2=sum(deep&np.isfinite(pos))
    X=np.c_[np.ones(p.sum()),n[p],d.E_f[p],np.log1p(d.UMI[p])]
    beta=np.linalg.lstsq(X,f[p],rcond=None)[0][1] if n1 and n2 and np.linalg.matrix_rank(X)==4 else np.nan
    reason='ok' if min(n1,n2)>=20 else 'empty_mask' if not np.any(inside) else 'no_observed_interface' if not np.isfinite(dist).any() else 'no_interior_scores_with_distance' if not sum(p) else 'no_deep_interior' if not n2 else 'no_near_points' if not n1 else 'below_20_per_band'
    out=dict(section_id=sec,boundary=b,variant=variant,width=w,n_B1=n1,n_B2=n2,mask_n=sum(inside),n_domain=sum(domain),n_dist=np.isfinite(dist).sum(),score=value,adjusted=beta,n_p95_B1=n_p95_B1,n_p95_B2=n_p95_B2,P95_eligible=min(n_p95_B1,n_p95_B2)>=20,P95_pp=prob if variant!='border_spot' else np.nan,eligible=min(n1,n2)>=20,reason=reason)
    effects.append(out);return out
for sec,r in M.iterrows():
    d=load(f'spots/{sec}.parquet');q=d.new_epi_domain.to_numpy();f=d.MP10.to_numpy();mu=avg(f[q])
    common=q&d[B[1]].to_numpy()&np.isfinite(d.dist_all_E)&np.isfinite(d['dist_'+B[1]])
    overall.append(dict(section_id=sec,rule='all',n=sum(q),mean_MP10=mu,median_MP10=np.median(f[q]) if any(q) else np.nan))
    for b in B:
        inside=d[b].to_numpy();dist=d['dist_'+b].to_numpy();bd=d['border_spot_'+b].to_numpy()
        tag=band(dist,200);ok=np.isfinite(dist);nt=band(bd,200)
        # A border center has zero distance on both sides. Retain explicit side.
        nt=np.where(ok&~inside&(abs(bd)<=200),'G1',nt)
        geos.append(dict(section_id=sec,boundary=b,mask_n=sum(inside),distance_n=sum(ok),inside_scored=sum(q&inside),max_abs_delta=np.max(abs(dist[ok]-bd[ok])) if any(ok) else np.nan,median_abs_delta=avg([]) if not any(ok) else np.median(abs(dist[ok]-bd[ok])),band_agreement=safe(sum(ok&(nt==tag)),sum(ok)),strict_n=d['strict_'+b].notna().sum()))
        for w in [100,200,300]:
            res=contrast(d,q,b,dist,w,'native',inside)
            for sl in ['B1','B2','G1','G2']:
                mask=band(dist,w)==sl
                for name,domain in [('E_scoring',q),('EIN_EIS',d.region_main.isin(['EIN','EIS'])),('non_E',~d.all_E&d.region_main.ne('Unknown'))]:
                    sel=mask&domain
                    counts.append(dict(section_id=sec,boundary=b,width=w,band=sl,domain=name,n=sum(sel),MP10=avg(f[sel]) if name!='non_E' else np.nan,epi_fraction=d.E_f[sel].mean(),T_fraction=d.T_f[sel].mean(),I_fraction=d.I_f[sel].mean(),UMI=d.UMI[sel].mean(),nGene=d.nGene[sel].mean(),top20_measurable=d.available_top20[sel].mean(),top20_detected=d.top20_detected[sel].mean()))
        for variant,dom,ds,mask in [('strict',q,d['strict_'+b].to_numpy(),inside),('border_spot',q,bd,inside),('legacy_domain',d.legacy_tumor_epi,dist,inside),('epi_query_subset',q&d.legacy_query,dist,inside),('common_points',common,dist,inside),('minus20',d.expression_QC&np.isfinite(f)&d.region_minus20.str.startswith('E'),d['minus_'+b].to_numpy(),d.region_minus20.str.startswith('E') if b==B[0] else d.region_minus20.isin(['ENN','ENS']))]:
            contrast(d,dom,b,ds,200,variant,mask)
        cc=[]
        for cid,a in d[q&inside&ok].groupby(b+'_component'):
            t=tag[a.index];v=a[t=='B1'];z=a[t=='B2'];wt=min(len(v),len(z))
            row=dict(section_id=sec,boundary=b,component=cid,B1=len(v),B2=len(z),score=v.MP10.mean()-z.MP10.mean(),weight=wt if wt>=3 else 0)
            comps.append(row);cc.append(row)
        cs=pd.DataFrame(cc,columns=['weight','B1','B2','score']);cs=cs[cs.weight>0];v=sum(cs.B1);z=sum(cs.B2)
        effects.append(dict(section_id=sec,boundary=b,width=200,variant='same_component',n_B1=v,n_B2=z,score=np.average(cs.score,weights=cs.weight) if len(cs) else np.nan,eligible=min(v,z)>=20,reason='ok' if min(v,z)>=20 else 'insufficient_common_component_support'))
        keep=q&ok&(tag!='B1')
        overall.append(dict(section_id=sec,rule='exclude_B1_'+b,n=sum(keep),mean_MP10=avg(f[keep]),median_MP10=np.median(f[keep]) if keep.sum()>=20 else np.nan))
        for x,a in d[q&ok].groupby(np.floor(dist[q&ok]/100)*100+50):
            bins.append(dict(section_id=sec,boundary=b,center=x,n=len(a),score=a.MP10.mean(),centered=a.MP10.mean()-mu,P95=a.P95_new.mean(),side='inside' if x<0 else 'EIN_EIS'))
        assert not (q&ok&(dist>0)).any() if b==B[0] else not (q&ok&(dist>0)&~d.region_main.isin(['EIN','EIS'])).any()
        if b==B[0]:both_native=inside
    for cat in CATS+['Unknown']:
        z=d.region_main.eq(cat);valid=d.region_main.ne('Unknown');inside=d.all_E
        areas.append(dict(section_id=sec,type='region',category=cat,n=sum(z),denominator=sum(valid) if cat!='Unknown' else len(d),fraction=safe(z.sum(),sum(valid) if cat!='Unknown' else len(d)),within_E=safe(z.sum(),sum(inside)) if cat.startswith('E') else np.nan,score=d.MP10[z&d.new_epi_domain].mean()))
    for k in range(9):
        z=d.accepted_label.eq(k);n=d.accepted_label.notna().sum();areas.append(dict(section_id=sec,type='R12',category=str(k),n=z.sum(),denominator=n,fraction=safe(z.sum(),n),score=d.MP10[z&d.new_epi_domain].mean()))
    for cl,a in d[d.cluster_main>=0].groupby('cluster_main'):
        clusters.append(dict(section_id=sec,local_cluster=cl,region=a.region_main.iloc[0],n=len(a),E_fraction=a.E_f.mean(),T_fraction=a.T_f.mean(),MP10_E=a.MP10[a.new_epi_domain].mean(),marker_file=str(V/'mia'/f'{sec}_markers.tsv.gz')))
e=pd.DataFrame(effects)
for variant in ['native','common_points']:
    for sec,d in e.query("width==200 and variant==@variant").groupby('section_id'):
        z=d.set_index('boundary');eligible=z.eligible.all();paired.append(dict(section_id=sec,domain=variant,eligible=eligible,D_all=z.at[B[0],'score'],D_without=z.at[B[1],'score'],paired_difference=z.at[B[1],'score']-z.at[B[0],'score'],reason_all=z.at[B[0],'reason'],reason_without=z.at[B[1],'reason']))
for key,rows in dict(section_effects=e,band_coverage=counts,border_distance_agreement=geos,component_contrasts=comps,section_bins=bins,section_exposures=overall,region_niche=areas,paired_boundaries=paired,within_section_clusters=clusters).items():save(rows,'tables/'+key+'.tsv.gz')

