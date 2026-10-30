#!/usr/bin/env python3
"""Post-hoc common exact-complex LR sensitivity from frozen full-background networks.
No inference, expression normalization, Q modification, or cache duplication occurs.
"""
from summary_common import *
from scipy import sparse
import argparse, hashlib, subprocess, traceback, time

METHODS=['COMMOT','SpatialCellChat'];STRATA=['dataset','tissue_condition','treatment_status','score_coverage_status']

def canonical(ligand,receptor,delimiter):
    def side(value):
        assert pd.notna(value), 'Missing resolved LR subunits'
        parts=[v.strip() for v in str(value).split(delimiter) if v.strip()]
        assert parts and all(';' not in v for v in parts)
        # Retain multiplicity and ligand/receptor direction; complexes are not split into single-gene pairs.
        return ';'.join(sorted(parts))
    return side(ligand)+'->'+side(receptor)

def lr_mapping(method,source):
    if method=='COMMOT':
        z=read(source/'position_lr_index.tsv');z['native_id']=z.lr_index.astype(str)
        z['canonical_LR']=[canonical(l,r,'_') for l,r in zip(z.ligand,z.receptor)]
        z['ligand_subunits']=[';'.join(sorted(str(v).split('_'))) for v in z.ligand]
        z['receptor_subunits']=[';'.join(sorted(str(v).split('_'))) for v in z.receptor]
    else:
        used=read(source/'position_LR_support.tsv')[['interaction']]
        cov=read(source/'gene_lr_coverage.tsv')
        assert not cov.interaction.duplicated().any(), 'Duplicate native SCC coverage interaction'
        z=used.merge(cov,on='interaction',how='left',validate='one_to_one');z['native_id']=z.interaction
        assert z.ligand_subunits.notna().all() and z.receptor_subunits.notna().all(), 'Unresolved SCC complex'
        z['canonical_LR']=[canonical(l,r,';') for l,r in zip(z.ligand_subunits,z.receptor_subunits)]
    assert not z.native_id.duplicated().any()
    return z

def write_commot(job,selected,out):
    inp,source=job['inp'],job['source'];meta=job['meta'];support=job['support'];n=len(meta)
    cross=sparse.csr_matrix((n,n));diag=np.zeros(n)
    for row in selected.itertuples():
        mat=sparse.load_npz(source/row.file).tocsr();assert mat.shape==(n,n) and np.isfinite(mat.data).all() and (mat.data>=0).all()
        diag+=mat.diagonal();mat.setdiag(0);mat.eliminate_zeros();cross+=mat
    for level in ['major','fine']:
        qdf=read(inp/f'Q_{level}.tsv').set_index('spot_id');assert qdf.index.tolist()==meta.spot_id.tolist()
        Q=qdf.to_numpy();N=qdf.columns.tolist();F=Q.T@cross@Q;D=np.outer(Q.sum(0),Q.sum(0))-Q.T@Q
        orig=read(source/f'{level}_edges.tsv');oldD=orig.pivot(index='source',columns='target',values='denominator').reindex(index=N,columns=N).to_numpy()
        assert np.allclose(D,oldD,rtol=1e-10,atol=1e-8,equal_nan=False), 'Frozen full-pair D differs'
        # Use exactly the original complete-pair denominator, including distance-induced zero pairings.
        D=oldD;G=np.divide(F,D,out=np.full_like(F,np.nan),where=D>0);same=Q.T@(diag[:,None]*Q)
        H=(Q>0).astype(float);binary=cross.copy();binary.data[:]=1;K=H.T@binary@H
        assert np.isclose(F.sum(),cross.sum()) and np.isclose(same.sum(),diag.sum())
        sx=support[support.level==level].set_index('node');good=sx.type_evaluable.reindex(N).map(truth).to_numpy();mask=good[:,None]&good[None,:]&(D>0)&bool(len(selected))
        rows=[]
        for i,s in enumerate(N):
            for j,t in enumerate(N):
                ok=bool(mask[i,j]);rows.append(dict(source=s,target=t,strength=G[i,j] if ok else np.nan,total_flow=F[i,j] if ok else np.nan,denominator=D[i,j],same_spot_flow=same[i,j] if ok else np.nan,positive_cross_links=K[i,j],computed=ok,n_LR=len(selected)))
        save(rows,out/f'{level}_edges.tsv')
        for key,M in [('total_strength',G),('total_flow',F),('denominator',D)]:
            pd.DataFrame(np.where(mask,M,np.nan),index=N,columns=N).rename_axis('node').to_csv(out/f'{level}_{key}.tsv',sep='\t',na_rep='NA')

def make_contrasts(job,out,ncommon,nstratum,stratum_id):
    info,section,method=job['info'],job['section'],job['method'];support=job['support'];cw=job['crosswalk'];mapping=job['mapping'];source=job['source']
    labelmap=dict(zip(cw.dropna(subset=['analysis_fine_label']).analysis_fine_label,cw.dropna(subset=['analysis_fine_label']).major))
    base={k:info.get(k,np.nan) for k in ['dataset','patient_uid','patient_id_status','analysis_unit_uid','analysis_unit_type','tissue_condition','treatment_status','score_coverage_status','hl_status','n_spots']}
    base.update(section_id=section,method=method,n_LR=len(mapping),n_original_LR=len(mapping),n_common_LR=ncommon,n_original_canonical_LR=mapping.canonical_LR.nunique(),common_stratum=stratum_id,n_stratum_sections=nstratum,common_LR_fraction=ncommon/len(mapping) if len(mapping) else np.nan,source_cache=str(source),sensitivity='posthoc_common_exact_complex_LR_no_refit')
    rows=[]
    for level in ['major','fine']:
        e=read(out/f'{level}_edges.tsv').rename(columns={'total_mean_strength':'strength','evaluable':'computed'}).set_index(['source','target']);orig=read(source/f'{level}_edges.tsv').rename(columns={'total_mean_strength':'strength'}).set_index(['source','target'])
        sx=support[support.level==level].set_index('node');good=sx.type_evaluable.map(truth);hlok=all(bool(good.get(n,False)) for n in EPI[:2]);targets=[n for n in sx.index if n not in EPI]
        for direction in ['epi_to_TME','TME_to_epi']:
            for node in targets:
                vals={};original={}
                for state,enode in [('high',EPI[0]),('low',EPI[1])]:
                    edge=(enode,node) if direction=='epi_to_TME' else (node,enode);er=e.loc[edge];vals[state]=er
                    original[state]=orig.loc[edge]
                ok=bool(hlok and good.loc[node] and ncommon>0 and all(pd.notna(vals[k].strength) for k in ['high','low']))
                hi,lo=(float(vals['high'].strength),float(vals['low'].strength)) if ok else (np.nan,np.nan)
                zero='NE' if not ok else 'both_zero' if hi==lo==0 else 'low_zero_only' if lo==0 else 'high_zero_only' if hi==0 else 'both_positive'
                ratio=np.log2(hi/lo) if zero=='both_positive' else np.nan
                oh,ol=float(original['high'].strength),float(original['low'].strength);original_ratio=np.log2(oh/ol) if ok and oh>0 and ol>0 else np.nan
                row={**base,'network_level':level,'lineage':node if level=='major' else labelmap.get(node,'unresolved'),'node':node,'direction':direction,'computed':ok,'NE_reason':'NA' if ok else 'empty_common_LR' if ncommon==0 else 'insufficient_state_or_type_support_or_denominator','strength_high':hi,'strength_low':lo,'delta_strength':hi-lo,'log2_HL':ratio,'zero_category':zero,'effect_direction':'NE' if not ok else 'higher' if hi>lo else 'lower' if hi<lo else 'equal','original_strength_high':oh if ok else np.nan,'original_strength_low':ol if ok else np.nan,'original_log2_HL':original_ratio,'delta_log2_HL_from_original':ratio-original_ratio,'direction_agrees_original':bool(np.sign(hi-lo)==np.sign(oh-ol)) if ok else pd.NA,'status':'computed_posthoc' if ok else 'NE'}
                for state in ['high','low']:
                    row[f'total_flow_{state}']=vals[state].total_flow if ok else np.nan
                rows.append(row)
    result=pd.DataFrame(rows)
    for native,common in [('strength_low','strength_low_common'),('strength_high','strength_high_common'),('log2_HL','log2_HL_common'),('n_common_LR','common_LR_n'),('n_original_LR','original_LR_n'),('common_LR_fraction','retained_fraction')]:result[common]=result[native]
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sections',nargs='+');p.add_argument('--output-root',type=Path,default=O/'comparison');p.add_argument('--mapping-only',action='store_true');p.add_argument('--force',action='store_true');a=p.parse_args();root=a.output_root;root.mkdir(parents=True,exist_ok=True)
    sections=a.sections or sorted(d.name for d in (O/'input_review/sections').iterdir() if (d/'READY.json').exists())
    jobs=[];failures=[];coverage=[]
    for section in sections:
        try:inp,info,meta,support,crosswalk=section_info(section)
        except Exception as exc:failures.append(dict(section_id=section,method='input',status='failed',reason=str(exc)));continue
        # Retain unknown tissue/treatment as distinct strata, never infer tumor from a name.
        info['score_coverage_status']='complete20' if int(info.get('n_mp10_genes_measured',0))==20 else 'score_coverage_incomplete'
        for method in METHODS:
            try:
                out,source,variant=job_dirs(method,section)
                if not ((source/'COMPLETED.json').exists() or (source/'COMPLETED.txt').exists()):
                    coverage.append(dict(section_id=section,method=method,status='not_complete',source_cache=str(source)));continue
                mapping=lr_mapping(method,source);jobs.append(dict(section=section,method=method,inp=inp,info=dict(info),meta=meta,support=support,crosswalk=crosswalk,source=source,mapping=mapping))
                coverage.append(dict(section_id=section,method=method,status='available',source_cache=str(source),n_native_LR=len(mapping),n_canonical_LR=mapping.canonical_LR.nunique()))
            except Exception as exc:failures.append(dict(section_id=section,method=method,status='failed',reason=str(exc),traceback=traceback.format_exc()))
    groups={}
    for job in jobs:groups.setdefault(tuple([job['method']]+[unknown(job['info'].get(k)) for k in STRATA]),[]).append(job)
    manifests=[];results=[];mapping_audits=[]
    for key,group in groups.items():
        common=set.intersection(*(set(j['mapping'].canonical_LR) for j in group));stratum_id=hashlib.sha256(json.dumps(key).encode()).hexdigest()[:16]
        for lr in sorted(common):manifests.append(dict(common_stratum=stratum_id,method=key[0],**dict(zip(STRATA,key[1:])),canonical_LR=lr,n_sections=len(group),section_ids=';'.join(j['section'] for j in group)))
        for job in group:
            section,method=job['section'],job['method'];mapping=job['mapping'];selected=mapping[mapping.canonical_LR.isin(common)].copy();out=root/'common_lr'/method/section;out.mkdir(parents=True,exist_ok=True)
            save(mapping,out/'canonical_LR_mapping.tsv');save(selected,out/'retained_native_LR.tsv')
            for other in jobs:
                if other['method']==method or other['section']!=section:continue
                inter=set(mapping.canonical_LR)&set(other['mapping'].canonical_LR)
                mapping_audits.append(dict(section_id=section,method=method,other_method=other['method'],n_canonical_LR=mapping.canonical_LR.nunique(),n_exact_complex_cross_method_overlap=len(inter)))
            identity=dict(section_id=section,method=method,stratum=list(key),sections=sorted(j['section'] for j in group),canonical_LR=sorted(common),selected_native_ids=selected.native_id.tolist(),source_cache=str(job['source']),source_completion_mtime=max((f.stat().st_mtime for f in [job['source']/'COMPLETED.json',job['source']/'COMPLETED.txt'] if f.exists()),default=0),support_mtime=(job['inp']/'type_support.tsv').stat().st_mtime,script_version=1)
            signature=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest();marker=out/'POSTHOC_COMPLETED.json'
            try:
                if a.mapping_only:continue
                reusable=marker.exists() and json.loads(marker.read_text()).get('signature')==signature and (out/'contrasts.tsv').exists()
                if reusable and not a.force:results.append(read(out/'contrasts.tsv'));print('COMMON LR REUSE',method,section,flush=True);continue
                if marker.exists():marker.unlink()
                if method=='COMMOT':write_commot(job,selected,out)
                else:
                    cmd=['/usr/local/bin/Rscript',str(O/'scripts/scc_common_lr.R'),section,str(out/'retained_native_LR.tsv'),str(out),'--source',str(job['source'])]
                    with (out/'aggregation.log').open('w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
                table=make_contrasts(job,out,len(selected),len(group),stratum_id);save(table,out/'contrasts.tsv');results.append(table)
                atom_json(dict(signature=signature,status='complete' if len(selected) else 'NE_empty_common_LR',n_native_LR=len(mapping),n_common_native_LR=len(selected),n_common_canonical_LR=len(common),n_stratum_sections=len(group),original_competition_background=True,no_refit=True,definition='COMMOT: sum retained native full-background LR then divide by original D; SCC: sum retained per-LR F/K',**identity),marker)
                print('COMMON LR',method,section,len(selected),'/',len(mapping),'stratum_n',len(group),flush=True)
            except Exception as exc:failures.append(dict(section_id=section,method=method,status='failed',reason=str(exc),traceback=traceback.format_exc(),recovery='rerun scripts/common_lr.py; completed matching strata reused'))
    save(coverage,root/'common_LR_input_coverage.tsv');save(pd.DataFrame(manifests,columns=['common_stratum','method',*STRATA,'canonical_LR','n_sections','section_ids']),root/'common_LR_retained_manifest.tsv');save(mapping_audits,root/'common_LR_cross_method_mapping_audit.tsv');save(pd.DataFrame(failures,columns=['section_id','method','status','reason','traceback','recovery']),root/'common_LR_failures.tsv')
    if not a.mapping_only:save(pd.concat(results,ignore_index=True) if results else pd.DataFrame(columns=['section_id','method','network_level','node','direction','computed']),root/'common_LR_sensitivity.tsv')
    atom_json(dict(status='complete_with_failures' if failures else 'mapping_only' if a.mapping_only else 'complete_available',n_requested_sections=len(sections),n_completed_native_jobs=len(jobs),n_strata=len(groups),n_failures=len(failures),n_not_complete=sum(r['status']=='not_complete' for r in coverage),native_inference_refitted=False,stratum_fields=STRATA,empty_intersection='NE; no artificial zero',resume='rerun exact command; changed cohort or retained set recomputes small aggregate only'),root/'COMMON_LR_STATUS.json')
    if failures:print('Recorded',len(failures),'failures; other jobs continued',flush=True)

if __name__=='__main__':main()
