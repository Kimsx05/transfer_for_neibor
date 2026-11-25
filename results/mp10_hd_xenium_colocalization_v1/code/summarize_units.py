from base import *
from display import *
from scipy.stats import spearmanr,pearsonr
from itertools import combinations
tax=pd.read_csv(O/'shared/taxonomy.tsv',sep='\t').drop_duplicates('label').set_index('label')
major=json.loads((O/'shared/major_markers.json').read_text())
sets=json.loads((O/'shared/programs.json').read_text())
focus={'T/NK','Endothelial','Fibroblast','Myeloid'}
def count_genes(x,genes,group):
    cols=pd.Index(genes).get_indexer(group); cols=cols[cols>=0]
    a=x[:,cols];return a.getnnz(axis=1),np.asarray(a.sum(1)).ravel()
def frozen(s):
    d=O/'hd'/s;e=pd.read_parquet(F/s/'Epi64_regions.parquet',columns=['region_id','x_um','y_um','measurement_eligible','MP10','H05'])
    spatial(s,e,'Epi64_frozen_MP10',e.MP10.where(e.measurement_eligible),'Frozen MP10',(0,1))
    colors.update({'H05':'#aa1645','other Epi':'#e4b286','ineligible':'#c9c9c9'})
    spatial(s,e,'Epi64_H05',np.where(e.measurement_eligible,np.where(e.H05,'H05','other Epi'),'ineligible'))
def hd(s):
    d=O/'hd'/s;z=pd.read_csv(d/'units.tsv.gz',sep='\t').set_index('region_id');x,genes,ids=read_export(d/'counts');assert np.array_equal(z.index,ids)
    fs=sorted((d/'weights').glob('[0-9][0-9][0-9][0-9][0-9].tsv.gz'));w=pd.concat([pd.read_csv(p,sep='\t').set_index('region_id') for p in fs]);assert len(w)==int((d/'fit_complete.txt').read_text()) and w.index.is_unique
    w=w.rename(columns={k.replace('/','_'):k for k in tax.index});assert set(w)==set(tax.index)
    w=w.div(w.sum(axis=1),axis=0);w.to_parquet(d/'weights.parquet');mw=w.T.groupby(tax.major).sum().T
    for a,keys in [(w,['candidate','support','second','runner']),(mw,['major_candidate','major_support','major_second','major_runner'])]:
        order=np.argsort(-a.values,axis=1);v=pd.DataFrame({keys[0]:a.columns.values[order[:,0]],keys[1]:a.values[np.arange(len(a)),order[:,0]],keys[2]:a.columns.values[order[:,1]],keys[3]:a.values[np.arange(len(a)),order[:,1]]},index=a.index);z=z.join(v)
    z['major']='Unknown';z['flag']='uncertain';z['supported']='Unknown';z['marker_detected']=np.nan;z['anchor_detected']=np.nan
    for k,markers in major.items():
        j=np.where(z.major_candidate.eq(k).fillna(False))[0];n,_=count_genes(x[j],genes,markers);u=z.iloc[j]
        mix=(u.major_runner>=.25)&((u.major_support-u.major_runner)<.25);good=(n>=1)&(u.major_support>=.50)&((u.major_support-u.major_runner)>=.12)&~mix
        z.loc[u.index[mix],['major','supported','flag']]=['Mixed','Mixed','mixed'];z.loc[u.index[good],['major','supported','flag']]=[k,k,'major_only']
    for lab,row in tax.iterrows():
        if row.major not in focus:continue
        j=np.where((z.major==row.major)&z.candidate.eq(lab).fillna(False))[0];u=z.iloc[j];n,_=count_genes(x[j],genes,row.markers.split());anc,_=count_genes(x[j],genes,row.anchors.split())
        z.loc[u.index,'marker_detected']=n;z.loc[u.index,'anchor_detected']=anc
        good=(n>=1)&(anc>=1)&(u.support>=.4)&((u.support-u.runner)>=.12)
        if lab=='Treg':good &= (n>=2)|(count_genes(x[j],genes,['FOXP3'])[0]>0)
        z.loc[u.index[good],['supported','flag']]=[lab,'subtype_supported']
    z.index.name='id';return z,x,genes
def xen(s):
    d=O/'xenium'/s;z=pd.read_parquet(d/'cells.parquet').rename(columns={'x_centroid':'x_um','y_centroid':'y_um'});z.index.name='id'
    for f in ['identity','state_transfer','ucell','expression_clusters']:
        t=pd.read_csv(d/(f+'.tsv.gz'),sep='\t').set_index('id');assert t.index.is_unique;z=z.join(t,validate='one_to_one')
    p=next((P/'highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5').glob('*_'+s+'_*h5'));x,ids,g=tenx(p); assert np.array_equal(ids,z.index);j=(g.feature_type=='Gene Expression').to_numpy();x=x[:,j];genes=g.name[j].to_numpy()
    ge=pd.Index(genes);panel=sorted(set(sum(sets.values(),[])));gg=ge.get_indexer(panel);assert min(gg)>=0
    raw=pd.DataFrame(x[:,gg].toarray(),index=ids,columns=panel); eligible=z.epi_eligible.fillna(False).to_numpy();ei=z.index[eligible]
    assert np.array_equal(z.Chol_biosynthesis_8.notna(),eligible) and np.array_equal(z.H05_transfer.notna(),eligible)
    diag=[]
    for k,ls in sets.items():
        a=raw.loc[ei,ls];z.loc[ei,k+'_detected']=(a>0).sum(1);z.loc[ei,k+'_max_gene_fraction']=a.max(axis=1)/a.sum(1).replace(0,np.nan)
        for gene in ls:diag.append((k,gene,len(a),int((a[gene]>0).sum()),int(a[gene].sum()),int(((a[gene]>0)&((a>0).sum(1)==1)).sum())))
    put(pd.DataFrame(diag,columns=['program','gene','epi','detected_epi','transcripts','sole_detected_gene_epi']),d/'epi_detection.tsv')
    v=z.loc[ei,'Chol_biosynthesis_8'];q=v.quantile(.95);high=v>q
    valid_group=bool(high.sum()>=5 and v.nunique()>2 and v.max()>0 and (z.loc[ei[high],'Chol_biosynthesis_8_detected']>=2).mean()>=.5)
    z['chol_group']=pd.NA
    if valid_group:z.loc[ei,'chol_group']=np.where(high,'Chol-high Epi','other Epi')
    put(pd.DataFrame([dict(n=len(ei),q95=q,n_at_q95=(v==q).sum(),strict_high=high.sum(),n_distinct=v.nunique(),group_evaluable=valid_group,rank_ties='linear quantile; strict >')]),d/'threshold.tsv')
    keys=['Chol_biosynthesis_8','MP10_covered','H05_transfer','REACTOME_CHOLESTEROL_BIOSYNTHESIS','HALLMARK_CHOLESTEROL_HOMEOSTASIS']
    corr=[]
    for a,b in list(combinations(keys,2))+[(k,h) for k in keys for h in ['nGeneTranscripts','nFeature','cell_area']]:
        vs=z.loc[ei,[a,b]].dropna();corr.append((a,b,len(vs),vs[a].corr(vs[b]),vs[a].corr(vs[b],method='spearman')))
    put(pd.DataFrame(corr,columns=['a','b','n','Pearson','Spearman']),d/'correlations.tsv')
    for k in keys:spatial(s,z,'state_'+k,z[k],k,(0,1))
    for gene in ['SREBF1','SREBF2']:spatial(s,z,'raw_'+gene,raw[gene].where(z.epi_eligible),'Raw transcripts')
    scatter(s,z.loc[ei,'Chol_biosynthesis_8'],z.loc[ei,'H05_transfer'],'Chol biosynthesis (UCell)','H05 transfer support','chol_H05')
    a=z.loc[ei,keys].rank(pct=True);discord=(a.Chol_biosynthesis_8-a.H05_transfer).abs()>=.7
    z.loc[ei[discord],['x_um','y_um']+keys].join(raw).to_parquet(d/'discordant_expression.parquet')
    ct=z.groupby(['cluster','supported'],observed=True).size().rename('cells').reset_index();put(ct,d/'cluster_identity.tsv')
    for col in ['area_outlier','boundary_degenerate','flag']:
        pass
    z.loc[z.area_outlier|z.boundary_degenerate|z.flag.eq('mixed'),['x_um','y_um','cell_area','segmentation_method','boundary_vertices','flag','major_candidate','major_second']].to_parquet(d/'segmentation_review.parquet')
    return z,x,genes

def deliver(s):
    ishd=s in HD;d=O/('hd' if ishd else 'xenium')/s
    z,x,genes=hd(s) if ishd else xen(s);z.to_parquet(d/'labels.parquet')
    counts=z.groupby(['major','supported','flag'],dropna=False).size().rename('units').reset_index();put(counts,d/'label_counts.tsv')
    for k in ['major','supported','flag']:spatial(s,z,k,z[k].replace({'mixed':'Mixed','uncertain':'Unknown'}))
    spatial(s,z,'candidate_support',z.support,'Leading subtype weight' if ishd else 'Within-lineage transfer support',(0,1))
    for group,gs in {'Treg':['FOXP3','IL2RA'],'T_NK':['CD3D','NKG7'],'CAF':['COL11A1','DPT'],'Endothelial':['PECAM1','VWF'],'Myeloid':['C1QC','SPP1']}.items():
        covered=[g for g in gs if g in genes]
        if not covered: continue
        n,total=count_genes(x,genes,covered);spatial(s,z,'marker_'+group,np.log1p(total),'log1p raw '+' + '.join(covered))
if __name__=='__main__':
    for s in sys.argv[1:]:
        if s=='frozen':
            for h in HD:frozen(h)
        else:deliver(s)

