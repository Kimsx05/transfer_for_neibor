from common import *
from collections import defaultdict
T=Path('tables')
root=Path('cache/oracle')
panel=rd(T/'panel_members.tsv').query('P_clean').gene
program=rd(T/'program_members.tsv')
methods=[]; audits=[];qc=[];allres=[];programs=[];members=[];edges=[];paths=[];patients=[];stateindex=[]
base=rd('cache/base_edges.tsv.gz')
edges.append(base.assign(record_type='prior',analysis='hg38',state='all'))
def ranked(x):
 x=x.copy();ok=x.status.eq('EVALUABLE') & ~x.is_KO
 x['rank']=np.nan;x['rank_denominator']=int(ok.sum())
 x.loc[ok,'rank']=x.loc[ok,'absolute_delta'].rank(ascending=False,method='average')
 x['rank_fraction']=(x['rank']-1)/max(1,int(ok.sum())-1)
 return x
def patient_summary(p):
 ids=['ko','graph','steps','gene']
 if p.empty:return pd.DataFrame(columns=ids+['signed_delta','absolute_delta','sd_delta','n_patients','n_cells'])
 q=p.groupby(ids+['patient_id'],sort=False)[['cells','signed_sum','abs_sum','std_sum']].sum(min_count=1).reset_index()
 for nm,c in [('signed_delta','signed_sum'),('absolute_delta','abs_sum'),('sd_delta','std_sum')]:q[nm]=q[c]/q.cells
 return q.groupby(ids,as_index=False).agg(signed_delta=('signed_delta','median'),absolute_delta=('absolute_delta','median'),sd_delta=('sd_delta','median'),n_patients=('patient_id','nunique'),n_cells=('cells','sum'))
def add_ranks(df):
 return pd.concat([ranked(g) for _,g in df.groupby(['ko','graph','steps'],sort=False)],ignore_index=True)
def program_table(tab,job):
 for (ko,graph,step),g in tab.groupby(['ko','graph','steps']):
  g=g.set_index('gene')
  for name in program.program.unique():
   for mask in ['Full','M20']:
    ids=program[(program.program==name)&(program['mask']==mask)].gene
    z=g.reindex(ids);present=z.n_states_expression.fillna(0)>0
    ok=z[z.status.eq('EVALUABLE')&z.is_KO.eq(False)]
    valid=len(ok)>=3 and 2*len(ok)>=len(z)
    programs.append(dict(analysis=job,ko=ko,graph=graph,steps=step,program=name,mask=mask,n_frozen=len(ids),n_in_panel=present.sum(),n_evaluable=len(ok),status='EVALUABLE' if valid else 'NE',reason='' if valid else 'empty mask or coverage below gate',negative_fraction=(ok.signed_delta<0).mean() if valid else np.nan,positive_fraction=(ok.signed_delta>0).mean() if valid else np.nan,median_signed_delta=ok.signed_delta.median() if valid else np.nan,median_absolute_delta=ok.absolute_delta.median() if valid else np.nan,median_rank_fraction=ok.rank_fraction.median() if valid else np.nan))
    if ko in K[:2]:
     members.append(z.reset_index().assign(analysis=job,ko=ko,graph=graph,steps=step,program=name,mask=mask))
for job in sorted(p for p in root.iterdir() if p.name!='pilot'):
 frames=[];pp=[]
 for st in sorted(p for p in job.iterdir() if p.is_dir()):
  marker=st/'complete.json';info=json.loads(marker.read_text()) if marker.exists() else dict(status='failed',reason='see worker log')
  methods.append(dict(analysis=job.name,state=st.name,**{k:info.get(k) for k in ['cells','patients','status','reason']}))
  if info['status']!='complete':continue
  a=rd(st/'all_genes.tsv.gz');p=rd(st/'patient.tsv.gz');aa=rd(st/'nodes.tsv.gz').assign(analysis=job.name,state=st.name);audits.append(aa)
  qc.append(rd(st/'simulation_qc.tsv.gz').assign(analysis=job.name,state=st.name));ed=rd(st/'edges.tsv.gz').assign(record_type='fitted',analysis=job.name,state=st.name);edges.append(ed)
  frames.append(a.assign(state=st.name));pp.append(p)
  sr=add_ranks(a.merge(patient_summary(p),how='left',on=['ko','graph','steps','gene']))
  sr['p.value']=np.nan;sr['p.adj']=np.nan
  sf=st/'ranked.tsv.gz';save(sr,sf);stateindex.append(dict(analysis=job.name,state=st.name,path=str(sf),rows=len(sr),sha256=sha(sf)))
  paths.append(a[a.ko.isin(K[:2]) & a.gene.isin(E+K+list(program.gene.unique()))].rename(columns={'ko':'source','gene':'target'}).assign(record_type='path',analysis=job.name,state=st.name)[['source','target','record_type','analysis','state','graph','steps','path_length','reachable','status','reason']])
 if not frames:continue
 raw=pd.concat(frames);p=pd.concat(pp);ids=['ko','graph','steps','gene']
 def merge_reasons(s):return '|'.join(sorted(set(s.dropna().astype(str))))
 raw['ok']=raw.status.eq('EVALUABLE');cover=raw.groupby(ids,as_index=False).agg(n_states_expression=('state','size'),n_states_evaluable=('ok','sum'));reason=raw[ids+['reason']].dropna().drop_duplicates().groupby(ids).reason.agg(merge_reasons);cover=cover.join(reason,on=ids)
 wh=raw[['ko','graph','steps']].drop_duplicates().merge(pd.DataFrame(dict(gene=panel)),how='cross')
 out=wh.merge(cover,on=ids,how='left').merge(patient_summary(p),on=ids,how='left');out['is_KO']=out.gene==out.ko;out['status']=np.where(out.n_patients.fillna(0)>0,'EVALUABLE','NE');out['reason']=out.reason.fillna('not mapped/expressed in any eligible cluster');out.loc[out.status=='EVALUABLE','reason']=''
 out=add_ranks(out);out['analysis']=job.name;out['p.value']=np.nan;out['p.adj']=np.nan;allres.append(out);program_table(out,job.name)
 if job.name=='main':save(p,T/'celloracle_patient_cluster_sums.tsv.gz')
save(pd.DataFrame(methods),T/'celloracle_job_status.tsv');save(pd.DataFrame(stateindex),T/'celloracle_cluster_index.tsv');save(pd.concat(audits),T/'celloracle_node_audit.tsv');save(pd.concat(qc),T/'celloracle_simulation_qc.tsv')
save(pd.concat(edges+paths,ignore_index=True),T/'celloracle_edges_and_paths.tsv')
overview=pd.concat(allres,ignore_index=True);save(overview,T/'celloracle_all_genes.tsv.gz')
save(overview.query('rank<=100'),T/'celloracle_top100.tsv.gz');save(overview.query('rank<=50'),T/'celloracle_top50.tsv.gz');save(overview.query('rank<=20'),T/'celloracle_top20.tsv.gz')
save(pd.DataFrame(programs),T/'celloracle_program_response.tsv');save(pd.concat(members),T/'celloracle_program_members.tsv.gz')
cross=[]
for j in overview.analysis.unique():
 sj='clean' if j=='main' else j
 s=rd(Path('cache/jobs')/sj/'heat.tsv.gz');s=s[s.ko.isin(K)&~s.is_KO][['ko','gene','distance','directionScore','depth_score','pc1_score','rank','rank_fraction','p.adj','reachable']].rename(columns={'rank':'sc_rank','rank_fraction':'sc_rank_fraction','p.adj':'sc_q','reachable':'sc_reachable'})
 co=overview[(overview.analysis==j)&(overview.graph=='native')&(overview.steps==3)&overview.ko.isin(K[:2])&~overview.is_KO]
 u=s.merge(co,how='outer',on=['ko','gene']);u['analysis']=j
 unsupported=u.ko.isin(K[2:]);u.loc[unsupported,'reason']='non-TF source outside standard CellOracle';u['status']=u.status.fillna('NE');u['co_direction']=np.sign(u.signed_delta)
 common=u.distance.notna() & u.signed_delta.notna()
 u['sc_common_rank']=u[common].groupby('ko').distance.rank(ascending=False);u['co_common_rank']=u[common].groupby('ko').absolute_delta.rank(ascending=False);u['n_common']=u.ko.map(u[common].groupby('ko').size()).fillna(0)
 u['comparison']=np.select([~common,u.co_direction==0,np.sign(u.directionScore)==u.co_direction],['NE','zero','agree'],default='disagree');cross.append(u)
save(pd.concat(cross,ignore_index=True),T/'cross_method_summary.tsv')

