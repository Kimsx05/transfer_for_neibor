from lib import *
import ast,h5py,traceback
from scipy.stats import rankdata
src=ROOT/'MP10_spatial_multisection_top20_vs_core19_v1/code/run_benchmark.py'
ns=dict(np=np,rankdata=rankdata);n=next(n for n in ast.parse(src.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='score')
exec(compile(ast.Module(body=[n],type_ignores=[]),str(src),'exec'),ns)
top=read(V/'audit/top20.tsv').gene.tolist()
mapping=read(V/'audit/reference_mapping.tsv').set_index('reference').compartment
groups={'E':['Epithelial'],'I':mapping[mapping=='I'].index.tolist(),'T':[k for k in mapping.index if k.startswith(('CD4','CD8','Treg','Tfh'))],'CD8':[k for k in mapping.index if k.startswith('CD8')],'Treg':['Treg'],'myeloid':[k for k in mapping.index if k.endswith((' TAM',' Mono',' DC')) or k in ['cDC1','cDC2','pDC','Neutrophil']]}
def score_and_qc(d,r):
    with h5py.File(r.input_path) as f:
        def strings(k):
            v=f[k];return (v if isinstance(v,h5py.Dataset) else v['values']).asstr()[:]
        a=f['X'];x=sparse.csr_matrix((a['data'][:],a['indices'][:],a['indptr'][:]),shape=a.attrs['shape'])
        pos=pd.Index(strings('obs/_index')).get_indexer(d.barcode);assert min(pos)>=0
        x=x[pos];genes=pd.Index(strings('var/gene_symbol'))
    if not genes.is_unique:
        j,genes=pd.factorize(genes);x=(x@sparse.csr_matrix((np.ones(len(j)),(np.arange(len(j)),j)))).tocsr()
    x.eliminate_zeros();d['nGene']=x.getnnz(axis=1);d['UMI']=np.asarray(x.sum(1)).ravel()
    assert np.allclose(d.UMI,d.nCount_RNA)
    d['expression_QC']=d.center_qc & (d.nGene>=150)&(d.UMI>0)
    d['available_top20']=genes.isin(top).sum()
    assert np.array_equal(x[:,genes.isin(top)].getnnz(axis=1),d.top20_detected)
    need=(d.all_E&d.expression_QC)|d.legacy_query|(d.region_minus20.str.startswith('E')&d.expression_QC)
    score=d[SCORE].copy();pre=np.isfinite(score);source=np.where(pre,'v1_stored','not_requested').astype(object)
    missing=need&~pre
    if missing.any():
        pp=ROOT/'neighborhood_v1/00_inputs_and_audit/mp_completion/by_section'/f'{r.section_id}_MP_Top20_AUCell.csv.gz'
        if pp.exists():
            older=pd.read_csv(pp,usecols=['barcode',SCORE]).set_index('barcode')[SCORE]
            revived=d.barcode.map(older);use=missing&revived.notna()
            score.loc[use]=revived[use];source[use]='upstream_stored'
    check=np.flatnonzero(d.legacy_query&pre)[::max(1,int((d.legacy_query&pre).sum()/12))][:12]
    if not len(check):check=np.flatnonzero(need&pre)[:12]
    def calc(idx):return ns['score'](x[idx],genes,top,int(np.ceil(.05*len(genes))))[0]
    err=max(abs(calc(check)-score.to_numpy()[check]),default=np.nan)
    assert not np.isfinite(err) or err<1e-7
    missing=need&score.isna()
    if missing.any() and d.available_top20.iloc[0]>0:
        score.loc[missing]=calc(np.flatnonzero(missing));source[missing]='recalculated'
    d['MP10']=score;d['score_source']=source;d['score_before']=pre;d['score_reason']=np.where(score.notna(),'available',np.where(~need,'not_requested','no_computable_signature'))
    d['new_epi_domain']=d.all_E&d.expression_QC&np.isfinite(score)
    target=score[d.new_epi_domain];cut=target.quantile(.95)
    d['P95_new']=pd.array(np.where(d.new_epi_domain,score>=cut,np.nan),dtype='Int8')
    return dict(n_check=len(check),max_error=err,n_genes=len(genes),rank_limit=int(np.ceil(.05*len(genes))),P95_threshold=cut,P95_fraction=d.P95_new.mean(),degenerate=bool(target.nunique()<=1 or cut==0),tie_excess=bool(d.P95_new.mean()>.10) if len(target) else False)
def one(sec,r):
    d=pd.read_parquet(V/'spots'/f'{sec}.parquet',use_threads=False)
    d['legacy_query']=d['query'];d['legacy_inside']=d.lesion_location=='inside';d['legacy_tumor_epi']=d.tumor_epi
    d['P95_old']=pd.array(np.where(d.tumor_epi,d.MP10_P95,np.nan),dtype='Int8')
    d['all_E']=d.region_main.str.startswith('E');d[B[1]]=d.region_main.isin(['ENN','ENS']);d['block']=d.tissue_component_id.astype('str')
    assert (~d[B[1]]|d.all_E).all()
    for w in [500,1000]:d['window'+str(w)]=d.block+'|'+(d.x_um//w).astype(int).astype(str)+'|'+(d.y_um//w).astype(int).astype(str)
    record=score_and_qc(d,r)
    c=pd.read_csv(r.c2l_mean_path,index_col=0).reindex(d.barcode)
    assert set(c)==set(mapping.index) and c.notna().all().all()
    d['total_abundance']=c.sum(axis=1).to_numpy()
    for k,members in groups.items():
        d[k+'_ab']=c[members].sum(axis=1).to_numpy();d[k+'_f']=d[k+'_ab']/d.total_abundance
    assert np.allclose(d.E_f,d.Epithelial_fraction,atol=1e-6)
    niche=pd.read_parquet(ROOT/'mp10_paper_distance_subtypes_v2/spots'/f'{sec}.parquet',columns=['barcode','accepted_label'],use_threads=False)
    d=d.merge(niche,on='barcode',validate='1:1')
    a=tissue_graph(d,d.region_main.ne('Unknown').to_numpy())
    coo=sparse.triu(a,k=1).tocoo();ij=coo.row,coo.col
    segs=[]
    for k,v1 in zip(B,['niche_all_epi','niche_low_immune']):
        mask=d[k].to_numpy();same=mask[ij[0]]==mask[ij[1]];pairs=(ij[0][same],ij[1][same])
        _,comp=connected_components(sparse.coo_matrix((np.ones(sum(same)),pairs),shape=a.shape),directed=False)
        d[k+'_component']=np.where(mask,comp,-1)
        seg=pd.read_parquet(V/'segments'/f'{sec}.parquet')
        if r.source_context=='adjacent_normal_source':
            fields,seg=geometry(d,d.region_main.astype(str).to_numpy(),v1)
            d['dist_'+k]=fields['d'];d['strict_'+k]=fields['strict']
            fields2,_=geometry(d,d.region_minus20.astype(str).to_numpy(),v1);d['minus_'+k]=fields2['d']
        else:
            seg=seg[seg.definition=='main_'+v1].copy()
            d['dist_'+k]=d['main_'+v1];d['strict_'+k]=d['main_'+v1+'_strict'];d['minus_'+k]=d['minus20_'+v1]
        border=np.full(len(d),np.nan)
        if len(seg):
            reg=d.set_index('barcode').region_main
            seg['outer_region']=seg.outside_barcode.map(reg).astype(str);seg['inner_region']=seg.inside_barcode.map(reg).astype(str)
            seg['contact_type']=np.where(seg.outer_region.isin(['EIN','EIS']),'to_EIN_EIS','to_non_E')
            for label,sgn in [('inside',-1),('outside',1)]:
                ids=d.barcode.isin(seg[label+'_barcode'])
                for component,v in d[ids].groupby(comp[ids]):
                    use=(comp==component)&d.region_main.ne('Unknown');border[use]=sgn*cKDTree(v[['x_um','y_um']]).query(d.loc[use,['x_um','y_um']])[0]
            seg['boundary']=k;seg['section_id']=sec;segs.append(seg)
        d['border_spot_'+k]=border
        d['status_'+k]=np.select([d.region_main=='Unknown',not mask.any(),d['dist_'+k].isna(),d['strict_'+k].isna()],['unknown','empty_mask','no_observed_contact_in_component','observed_only'],default='strict')
        assert np.array_equal(np.isfinite(border),np.isfinite(d['dist_'+k]))
        if not mask.any():assert not np.isfinite(d['dist_'+k]).any()
    if (d[B[0]]==d[B[1]]).all():assert np.allclose(d['dist_'+B[0]],d['dist_'+B[1]],equal_nan=True)
    keep=['section_id','barcode','array_row','array_col','x','y','x_um','y_um','block','window500','window1000','region_main','region_minus20','cluster_main','accepted_label','UMI','nGene','top20_detected','available_top20','MP10','score_before','score_source','score_reason','expression_QC','legacy_query','legacy_inside','legacy_tumor_epi','new_epi_domain','P95_new','P95_old','all_E',B[1],'total_abundance']+[k+v for k in A for v in ['_ab','_f']]+[col for k in B for col in [k+'_component','dist_'+k,'strict_'+k,'minus_'+k,'border_spot_'+k,'status_'+k]]
    pq(d[keep],f'spots/{sec}.parquet')
    boundary=pd.concat(segs,ignore_index=True) if segs else pd.DataFrame(columns=['section_id','boundary','inner_region','outer_region','contact_type','inside_barcode','outside_barcode','x1','x2','y1','y2'])
    pq(boundary,f'boundary/{sec}.parquet')
    cross=d.assign(score_after=d.MP10.notna()).groupby(['all_E','legacy_query','legacy_inside','score_before','score_after','expression_QC'],observed=True).size().reset_index(name='n')
    save(cross.assign(section_id=sec),f'audit/{sec}_cross.tsv')
    # Subtype means on shared domains without copying full posterior.
    rows=[]
    for domain,ok in [('whole',np.ones(len(d),bool)),('all_E',d.all_E)]:
        for name in c:
            vals=c[name].to_numpy()[ok];den=d.total_abundance[ok].sum()
            rows.append(dict(section_id=sec,domain=domain,cell_type=name,ab_sum=vals.sum(),all30_sum=den,n=sum(ok),fraction=safe(vals.sum(),den)))
    save(rows,f'metrics/{sec}_subtypes.tsv')
    return dict(record,section_id=sec,n=len(d),E=int(d.all_E.sum()),eligible=int(d.new_epi_domain.sum()),previously_shown=int((d.new_epi_domain&d.legacy_query&d.score_before).sum()),hidden_existing=int((d.new_epi_domain&~d.legacy_query&d.score_before).sum()),rescored=int((d.score_source=='recalculated').sum()),restored_upstream=int((d.score_source=='upstream_stored').sum()),zero=int((d.new_epi_domain&d.MP10.eq(0)).sum()),old_query=int(d.legacy_query.sum()),old_tumor_epi=int(d.legacy_tumor_epi.sum()))
if __name__=='__main__':
    out=[]
    for sec,r in M.iterrows():
        done=O/'logs'/f'{sec}.json'
        if done.exists() and all(path.is_file() and path.stat().st_size>0 for path in [O/'spots'/f'{sec}.parquet',O/'boundary'/f'{sec}.parquet',O/'audit'/f'{sec}_cross.tsv',O/'metrics'/f'{sec}_subtypes.tsv']):out.append(json.loads(done.read_text()));continue
        rec=one(sec,r);done.write_text(json.dumps(rec,default=lambda x:x.item()));out.append(rec);print(sec,flush=True)
    save(out,'audit/score_coverage_audit.tsv')
    allc=pd.concat([read(p) for p in (O/'audit').glob('*_cross.tsv')]);save(allc,'audit/score_cross_section.tsv');save(allc.groupby([k for k in allc if k not in ['n','section_id']]).n.sum().reset_index(),'audit/score_cross_cohort.tsv')
    save(M,'audit/input_manifest.tsv')

