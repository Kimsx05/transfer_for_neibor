from lib import *
EPS=1.0
def build(sec):
    d=load(f'spots/{sec}.parquet')
    d['tissue_component_id']=d.block
    valid=d.expression_QC.to_numpy() & d.block.notna().to_numpy() & np.isfinite(d.T_f)
    adj=tissue_graph(d,valid)
    all_e=d.all_E.to_numpy();dist=d.dist_all_E.to_numpy();xy=d[['x_um','y_um']].to_numpy()
    metrics=[]
    def add(name,num,den,domain,scale,formula,unit='ratio',value=None,reason=None):
        val=safe(num,den) if value is None else value
        metrics.append(dict(section_id=sec,metric=name,value=val,numerator=num,denominator=den,domain=domain,scale=scale,formula=formula,unit=unit,status=reason or ('ok' if np.isfinite(val) else 'undefined_or_empty_domain')))
    domains={'whole':valid,'inside':valid&all_e}
    for w in [100,200,300]:
        domains['outer'+str(w)]=valid&(dist>0)&(dist<=w)
        domains['transition'+str(w)]=valid&(abs(dist)<=w)
    for domain,mask in domains.items():
        n=sum(mask);total=d.total_abundance[mask].sum()
        for k in A:
            amt=d[k+'_ab'][mask].sum()
            add(f'{domain}_{k}_abundance',amt,1,domain,'section','sum posterior mean','estimated cells',value=amt if n else np.nan)
            add(f'{domain}_{k}_fraction',amt,total,domain,'section','sum posterior / sum all30')
            add(f'{domain}_{k}_density',amt,n*AREA,domain,'section','sum posterior / observed hex area','estimated cells/mm²')
    for k in A[1:]:
        inside=d[k+'_ab'][all_e].sum();total=d[k+'_ab'][valid].sum()
        add('entry_share_'+k,inside,total,'all_E/whole','section','inside abundance / whole abundance')
        add('area_corrected_entry_'+k,safe(inside,total),safe(sum(all_e),sum(valid)),'all_E/whole','section','entry share / E area share')
        for w in [100,200,300]:
            out=domains['outer'+str(w)];di=safe(inside,sum(all_e)*AREA);do=safe(d[k+'_ab'][out].sum(),sum(out)*AREA)
            add(f'penetration_{k}_{w}',di+EPS,do+EPS,'all_E/outer',str(w),'log2((inside density + 1)/(outer density + 1))','log2 ratio',value=np.log2((di+EPS)/(do+EPS)),reason='low_outer_support' if 0<sum(out)<5 else None)
    for order in ['same','R1','R12']:
        g=adj if order=='R1' else (adj+adj@adj).tocsr() if order=='R12' else sparse.diags(valid.astype(float)).tocsr()
        if order=='R12':g.setdiag(0);g.eliminate_zeros();g.data[:]=1
        degree=np.asarray(g.sum(1)).ravel();use=degree>0
        W=sparse.diags(1/np.maximum(degree,1))@g
        for k in ['T','I','CD8','E']:
            x=d[k+'_f'].to_numpy();z=(x-np.mean(x[use]))/np.std(x[use]) if np.std(x[use])>0 else np.full(len(d),np.nan)
            if order!='same':
                lag=W@z
                add(f'Moran_{k}_{order}',np.sum(z[use]*lag[use]),sum(use),'observed_nonisolates',order,'mean z(x)*Wz(x)')
                if k=='T' and order=='R1':
                    d['local_Moran_T']=np.where(use,z*lag,np.nan);d['lag_T']=np.where(use,W@x,np.nan)
            if k=='E':continue
            e=d.E_f.to_numpy();ex=(e-e[use].mean())/e[use].std()
            for a,b,y,t in [('E',k,e,x),(k,'E',x,e)]:
                cross=np.sum(y[use]*(W@t)[use]);denom=y[use].sum()
                add(f'co_{a}_{b}_{order}',cross,sum(use),'observed_nonisolates',order,'mean fraction A * W fraction B')
                add(f'CLQ_{a}_{b}_{order}',safe(cross,denom),t[use].mean(),'observed_nonisolates',order,'A-weighted local B fraction / global B fraction')
            if order!='same':
                add(f'contact_{k}_E_{order}',np.sum(x[use]*(W@e)[use]),np.sum(x[use]*(W@x)[use]),'observed_nonisolates',order,'sum T * WE / sum T * WT (replace T by named lineage)')
                add(f'biv_E_{k}_{order}',np.sum(ex[use]*(W@z)[use]),sum(use),'observed_nonisolates',order,'mean z(E)*Wz(B)')
                add(f'biv_{k}_E_{order}',np.sum(z[use]*(W@ex)[use]),sum(use),'observed_nonisolates',order,'mean z(B)*Wz(E)')
                add(f'Lee_E_{k}_{order}',np.sum((W@ex)[use]*(W@z)[use]),sum(use),'observed_nonisolates',order,'Lee L row-normalized nonisolates; mean Wz(E)*Wz(B)')
    # All local numerators are additive across observed sites in one section.
    s=pd.DataFrame(index=d.index)
    s['obs']=valid.astype(int);s['area_mm2']=valid*AREA;s['E_n']=all_e.astype(int)
    q=d.new_epi_domain.to_numpy();per=domains['outer200']
    s['per_n']=per.astype(int);s['scored_n']=q.astype(int)
    s['MP10_sum']=d.MP10.where(q,0);s['P95_n']=d.P95_new.fillna(0).astype(int)
    for label,mask in [('inside',all_e),('outer',per),('whole',valid)]:
        for k in A:s[label+'_'+k]=d[k+'_ab'].where(mask,0)
        s[label+'_total']=d.total_abundance.where(mask,0)
    s['lag_E_n']=np.isfinite(d.lag_T)&all_e;s['lag_T_E_sum']=d.lag_T.where(all_e,0).fillna(0)
    for k in B:
        for label in ['B1','B2']:
            use=q&(band(d['dist_'+k],200)==label)
            s[k+'_'+label+'_n']=use.astype(int);s[k+'_'+label+'_sum']=d.MP10.where(use,0);s[k+'_'+label+'_P95']=d.P95_new.where(use,0).fillna(0)
    locals=[]
    for key in ['block','window500','window1000']+[k+'_component' for k in B]:
        keep=valid if not key.endswith('component') else d[key].ne(-1)&valid
        z=s[keep].groupby(d.loc[keep,key],observed=True).sum().reset_index().rename(columns={key:'local_id'})
        z['local_id']=z.local_id.astype(str);z['level']=key;z['section_id']=sec
        z['E_T_density']=z.inside_T/(z.E_n*AREA).replace(0,np.nan);z['out_T_density']=z.outer_T/(z.per_n*AREA).replace(0,np.nan)
        z['P_T']=np.log2((z.E_T_density+EPS)/(z.out_T_density+EPS))
        z['E_T_fraction']=z.inside_T/z.inside_total.replace(0,np.nan);z['T_neighbor_E']=z.lag_T_E_sum/z.lag_E_n.replace(0,np.nan)
        z['MP10_mean']=z.MP10_sum/z.scored_n.replace(0,np.nan)
        z['entry_support']=np.select([z.E_n<5,z.per_n<4],['limited_E','limited_outer'],default='supported')
        locals.append(z)
        if key in ['window500','window1000','block']:
            for col in ['E_T_density','E_T_fraction','P_T','T_neighbor_E']:
                ok=(z.E_n>=5)&(z.per_n>=4 if col=='P_T' else True)
                v=z.loc[ok,col].dropna()
                for stat in ['IQR','median','MAD']:
                    val=(v.quantile(.75)-v.quantile(.25) if stat=='IQR' else v.median() if stat=='median' else (v-v.median()).abs().median()) if len(v)>=2 else np.nan
                    add(f'{key}_{col}_{stat}',val,1,'local units',key,stat+' of supported local '+col,value=val)
    for k in B:
        vc=d.loc[d[k],k+'_component'].value_counts();ds=d['dist_'+k];dd=abs(ds[d[k]&ds.notna()])
        for name,val in [('area',d[k].sum()*AREA),('fraction',d[k].mean()),('components',len(vc)),('median_component_n',vc.median()),('max_depth_um',dd.max()),('median_depth_um',dd.median()),('B1_n',sum(q&(band(ds,200)=='B1'))),('B2_n',sum(q&(band(ds,200)=='B2')))]:add(k+'_'+name,val,1,'structure','section',name,value=val)
    for name,val in [('mixed_E_share',safe(d.region_main.isin(['EIN','EIS']).sum(),sum(all_e))),('blocks',d.block.nunique()),('observed_area',valid.sum()*AREA),('UMI_median',d.UMI[valid].median()),('nGene_median',d.nGene[valid].median()),('score_coverage',safe(sum(q),sum(all_e)))]:add(name,val,1,'structure/QC','section',name,value=val)
    pq(pd.concat(locals,ignore_index=True),f'local/{sec}.parquet')
    pq(d[['section_id','barcode','local_Moran_T','lag_T']],f'metrics/{sec}_local_moran.parquet')
    return pd.DataFrame(metrics)
if __name__=='__main__':
    rows=[]
    for sec in M.index:
        p=f'metrics/{sec}_candidate.parquet'
        if not all(path.is_file() and path.stat().st_size>0 for path in [O/p,O/'local'/f'{sec}.parquet',O/'metrics'/f'{sec}_local_moran.parquet']):pq(build(sec),p)
        rows.append(load(p))
    data=pd.concat(rows);save(data,'tables/section_metrics_all.tsv.gz');save(data[['metric','domain','scale','formula','unit']].drop_duplicates(),'audit/metric_dictionary.tsv')
    local=pd.concat([load(f'local/{sec}.parquet') for sec in M.index]);pq(local,'tables/window_component_metrics.parquet')
    # Physically separate immune-only inputs, preserving structural NA.
    immune=data[~data.metric.str.contains('score|MP10|P95',case=False)]
    save(immune.pivot(index='section_id',columns='metric',values='value').reset_index(),'classification/immune_input.tsv.gz')
    save(local.loc[local.level.isin(['window500','window1000']),['section_id','level','local_id','obs','E_n','per_n','E_T_density','out_T_density','P_T','E_T_fraction','T_neighbor_E','area_mm2','inside_T','outer_T']],'classification/immune_windows.tsv.gz')

