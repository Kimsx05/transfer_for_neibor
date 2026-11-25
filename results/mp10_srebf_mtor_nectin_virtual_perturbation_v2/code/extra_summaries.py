from common import *
T=Path('tables');p=rd(T/'patient_conditioned.tsv.gz');m=rd(T/'program_members.tsv')
z=p[p.status=='EXPLORATORY'].merge(m[['gene','program','mask']],on='gene')
idx=['patient_id','ko','stratum','program','mask']
u=z.groupby(idx).agg(n_cells=('n_cells','first'),n_genes=('gene','nunique'),median_distance=('distance','median')).reset_index()
u=u.join(m.groupby(['program','mask']).size().rename('n_frozen'),on=['program','mask']);u['status']=np.where((u.n_genes>=3)&(2*u.n_genes>=u.n_frozen),'EXPLORATORY','NE');u.loc[u.status=='NE','median_distance']=np.nan
save(u,T/'patient_conditioned_programs.tsv.gz');cols=['patient_id','ko','program','mask'];good=u[u.status=='EXPLORATORY'];v=good[good.stratum=='L50'][cols+['median_distance','n_cells']].rename(columns={'median_distance':'low','n_cells':'n_low'})
h=good[good.stratum.isin(['H05','H10','H50'])].merge(v,on=cols);h['high_minus_low']=h.median_distance-h.low;save(h,T/'patient_paired_programs.tsv.gz')
# Program member values across resampling can be recovered by a keyed join to the full matrices.
# Record target/path direction coverage and historical control effects on current shared WT.
hist=rd(T/'legacy_controls.tsv');gm=rd(T/'gene_mapping.tsv').set_index('feature_id').gene;hist['current_gene']=hist.gene.map(gm)
a=rd('cache/jobs/clean/heat.tsv.gz');qq=hist[['gene','current_gene','matched_to','label']].rename(columns={'gene':'v1_control'}).merge(a[a.gene.isin(E)],left_on='current_gene',right_on='ko',how='left');save(qq,T/'legacy_control_endpoints.tsv')
import platform
env=dict(python=sys.version,platform=platform.platform(),analysis_python=sys.executable,reference_paths=[str(V),str(prior)],oracle_venv='cache/co310',r_tarballs={p.name:sha(p) for p in (V/'cache').glob('scTenifold*.tar.gz')},wheel={p.name:sha(p) for p in (V/'cache').glob('celloracle*.whl')})
Path('config/environment.json').write_text(json.dumps(env,indent=2)+'\n')

