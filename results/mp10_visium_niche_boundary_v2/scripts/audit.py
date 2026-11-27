from lib import *
import importlib.metadata as im,re
from frozen_geometry import nearest
st=read(O/'audit/source_studies.tsv').set_index('dataset_id');bc=read(O/'audit/gse319536_table_s1.tsv').set_index('specimen');clinical=[];checks=[];gates=[];periphery=[]
ident=ROOT/'cell2location_final_validation_gate_v1/09_final_status/DOWNSTREAM_USAGE_TABLE.csv'
save(pd.read_csv(ident),'audit/c2l_reference_usage.tsv')
from base import groups
save([dict(axis=k,member=m) for k,vs in groups.items() for m in vs],'audit/c2l_axes.tsv')
reg=read(ROOT/'neighborhood_ring_harmony_v3/04_comparison/niche_dictionary.tsv');save(reg[reg.config_id.eq('SAMPLE_R12_LEIDEN')],'audit/R12_dictionary.tsv')
ef=read(O/'tables/section_effects.tsv.gz').query("variant=='native' and width==200")
for sec,r in M.iterrows():
    ds=r.dataset_id;k=sec.split('_')[-1];s=st.loc[ds]
    clin=dict(section_id=sec,dataset_id=ds,doi=s.doi,source_sample_id=r.source_sample_id,patient_uid=r.patient_uid,analysis_unit_uid=r.analysis_unit_uid,source_context=r.source_context,evidence=r.evidence_method,local_stage='unknown',local_pretreatment='unknown',local_location='unknown',local_subsequent_regimen='unknown',published_TLS='unknown',published_immune_context=s.published_immune_observation)
    if ds=='GSE319536':
        j=bc.loc[k];clin.update(local_stage=j.stage,local_pretreatment='no prior treatment',local_location='primary bladder, diagnostic TURBT',local_subsequent_regimen=j.subsequent_NAC,published_TLS=j.TLS,source_table='Table S1, BC identifier')
    if ds=='GSE171351':clin.update(local_pretreatment='treatment-naive',local_stage='MIBC, exact T unknown',local_location='primary bladder')
    if ds=='GSE238145':clin.update(local_stage='invasive; exact TNM unknown',local_location='FFPE primary bladder')
    if ds=='GSE285715':clin.update(local_location='adjacent normal bladder' if '_N' in sec else 'tumor bladder')
    if ds.startswith('NPJ'):clin.update(local_stage='pT1 high-grade (all seven)',local_location='bladder en bloc TURBT; subsite unknown')
    clinical.append(clin)
    sib=pd.read_csv(Path(r.c2l_mean_path).parent/'section_sibling_correlation.csv');v=sib[(sib.family=='CD8_NK')&(sib.metric=='mean')&((sib.subtype_A=='NK')|(sib.subtype_B=='NK'))]
    rho=v.spearman.abs().max();gates.append(dict(section_id=sec,n_pairs=len(v),max_abs_CD8_NK_rho=rho,low_collinearity_screen=len(v)==2 and rho<.9))
    d=load(f'spots/{sec}.parquet');e=load(f'boundary/{sec}.parquet');d['tissue_component_id']=d.block;keys=d.set_index('barcode')
    equiv=(d[B[0]]==d[B[1]]).all()
    assert d.barcode.is_unique and not (d[B[1]]&~d.all_E).any() and not d.P95_new[~d.new_epi_domain].notna().any()
    for b in B:
        z=e[e.boundary==b];mask=d[b].to_numpy();dist=d['dist_'+b].to_numpy();valid=d.region_main.ne('Unknown').to_numpy()
        assert (keys.loc[z.inside_barcode,b]==True).all() and (~keys.loc[z.outside_barcode,b]).all()
        assert (keys.loc[z.inside_barcode,'block'].to_numpy()==keys.loc[z.outside_barcode,'block'].to_numpy()).all()
        assert not z.outer_region.eq('Unknown').any() and not np.isfinite(dist[~valid]).any()
        if not mask.any():assert not np.isfinite(dist).any() and not len(z)
        eff=ef[(ef.section_id==sec)&(ef.boundary==b)].iloc[0]
        if not len(z):assert not eff.eligible and not np.isfinite(eff.score)
        if equiv:
            assert np.allclose(d.dist_all_E,d['dist_'+B[1]],equal_nan=True)
            other=ef[(ef.section_id==sec)&(ef.boundary==B[1])].iloc[0]
            assert np.isclose(eff.score,other.score,equal_nan=True)
        checks.append(dict(section_id=sec,boundary=b,mask_exists=mask.any(),observed_interface=len(z)>0,segments=len(z),E_E_contact=(z.contact_type=='to_EIN_EIS').sum(),E_nonE_contact=(z.contact_type=='to_non_E').sum(),equal_masks=equiv,finite_distances=np.isfinite(dist).sum(),B1_n=eff.n_B1,B2_n=eff.n_B2,eligible=eff.eligible,reason=eff.reason,validation='pass'))
        # Allocate each outside G1 spot to a single adjacent E component, resolving ties in frozen segment order.
        adj=tissue_graph(d,valid);t=sparse.triu(adj,k=1).tocoo();same=mask[t.row]==mask[t.col]
        _,comp=connected_components(sparse.coo_matrix((np.ones(sum(same)),(t.row[same],t.col[same])),shape=adj.shape),directed=False)
        assert np.array_equal(comp[mask],d.loc[mask,b+'_component'])
        assigned=np.full(len(d),-1)
        for cc,v in z.groupby('outside_component'):
            ix=np.flatnonzero((comp==cc)&(dist>0)&(dist<=200));lines=v[['x1','y1','x2','y2']].to_numpy().reshape(-1,2,2)
            if len(ix):_,i=nearest(d.loc[ix,['x_um','y_um']].to_numpy(),lines);assigned[ix]=v.inside_component.to_numpy()[i]
        assert sum(assigned>=0)==sum((dist>0)&(dist<=200))
        for cc,a in d[mask].groupby(b+'_component'):
            rr=d[assigned==cc];record=dict(section_id=sec,boundary=b,component=cc,inner_n=len(a),unique_outer_n=len(rr),inner_area=len(a)*AREA,outer_area=len(rr)*AREA,inner_sum_all30=a.total_abundance.sum(),outer_sum_all30=rr.total_abundance.sum())
            for k in A:
                vi=a[k+'_ab'].sum();vo=rr[k+'_ab'].sum();record.update({k+'_inner_sum':vi,k+'_outer_sum':vo,k+'_penetration':np.log2((safe(vi,len(a)*AREA)+1)/(safe(vo,len(rr)*AREA)+1))})
            periphery.append(record)
save(clinical,'audit/section_clinical.tsv');save(checks,'module_status.tsv');save(gates,'audit/CD8_identifiability.tsv');save(periphery,'tables/component_unique_periphery.tsv.gz')
a=read(O/'tables/section_metrics_all.tsv.gz')
g=pd.DataFrame(gates)
save(a[a.metric.str.contains('CD8')].merge(g).query('low_collinearity_screen'),'tables/CD8_low_collinearity_subset.tsv.gz')
md=a[['metric','unit','domain','scale','formula']].drop_duplicates().copy()
for i,r in md.iterrows():
    key=r.metric;unit=r.unit
    if key.endswith('_area') or key=='observed_area':unit='mm²'
    if 'depth_um' in key:unit='μm'
    if key.endswith(('_n','_components')) or key=='blocks':unit='count'
    if key=='UMI_median':unit='UMI'
    if key=='nGene_median':unit='detected genes'
    if 'E_T_density_' in key:unit='estimated cells/mm²'
    if '_P_T_' in key:unit='log2 ratio'
    md.loc[i,'unit']=unit;md.loc[i,'direction']='higher interior relative to periphery' if key.startswith('penetration_') else 'greater same-lineage spatial clustering' if key.startswith('Moran_') else 'greater observed relative to expected neighborhood abundance' if key.startswith('CLQ_') else 'greater variability across local units' if key.endswith(('IQR','MAD')) else 'larger value of named measure'
md['missing_rule']='empty domain, zero denominator, insufficient local replicates or constant standardized variable → NA; consult row status'
md['epsilon']='none';md.loc[md.metric.str.contains('penetration|_P_T_'),'epsilon']='1 estimated cell/mm² added to both densities'
save(md,'audit/metric_dictionary.tsv')
a=a.drop(columns='unit').merge(md[['metric','unit']],on='metric',validate='m:1');save(a,'tables/section_metrics_all.tsv.gz')
save([dict(path=str(f),sha256=sha(f)) for f in [V/'scripts/common.py',V/'scripts/geometry.py',ROOT/'MP10_spatial_multisection_top20_vs_core19_v1/code/run_benchmark.py',ident]+list((O/'audit').glob('v1_*'))],'audit/software_inputs.tsv')
save([dict(package=p,version=im.version(p)) for p in ['numpy','pandas','scipy','scikit-learn','shapely','h5py','pyarrow','matplotlib']],'audit/software.tsv')
print('validated sections:',len(M),'geometry rows:',len(checks),'CD8 lower cross-lineage correlation:',sum(g.low_collinearity_screen))

