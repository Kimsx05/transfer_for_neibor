exec(open(__file__.replace('04b_patient_summary.py','01_prepare.py')).read().split('# Archive contents')[0])
from scipy.stats import binomtest
x=pd.read_csv(B/'01_epi/sample_effects.tsv',sep='\t');valid=x[x.patient_status.isin(['verified','confirmed'])&x.patient.notna()&x.rho.notna()].copy();out=[]
for keys,g in valid.groupby(['subset','dataset','patient','signature','adjustment']):
 rho=g.rho.to_numpy();finite=np.abs(rho)<1
 out.append(dict(zip(['subset','dataset','patient','signature','adjustment'],keys),n_samples=len(g),n_cells=int(g.n_used.sum()),fisher_z=float(np.mean(np.arctanh(rho))) if finite.all() else np.nan,reason='OK' if finite.all() else 'exact_abs_rho_1_nonfinite_Fisher_z'))
pdout=pd.DataFrame(out,columns=['subset','dataset','patient','signature','adjustment','n_samples','n_cells','fisher_z','reason']);pdout['rho']=np.tanh(pdout.fisher_z.astype(float));save(pdout,'01_epi/patient_effects.tsv');summ=[];rng=np.random.default_rng(cfg['seed'])
for subset in x['subset'].unique():
 for ds in ['ALL']+sorted(x.dataset.unique()):
  for sig in x.signature.unique():
   for adj in x.adjustment.unique():
    z=pdout[(pdout['subset']==subset)&(pdout.signature==sig)&(pdout.adjustment==adj)]
    if ds!='ALL':z=z[z.dataset==ds]
    # Patients with multiple samples are already merged and remain whole bootstrap units.
    assert not z.patient.duplicated().any()
    v=z.fisher_z.to_numpy(float);v=v[np.isfinite(v)];n=len(v);r=dict(subset=subset,dataset=ds,signature=sig,adjustment=adj,n_patients=n,rho=np.nan,CI_low=np.nan,CI_high=np.nan,P=np.nan,FDR=np.nan,status='BLOCKED_no_verified_patient_mapping')
    if n>=2:
     bt=np.tanh(rng.choice(v,(10000,n),replace=True).mean(axis=1));nz=v[v!=0];r.update(rho=np.tanh(v.mean()),CI_low=np.quantile(bt,.025),CI_high=np.quantile(bt,.975),P=binomtest((nz>0).sum(),len(nz),.5).pvalue if len(nz) else np.nan,status='OK')
    summ.append(r)
s=pd.DataFrame(summ)
# Only ALL/tumor original UP+DN unadjusted tests form the prespecified primary family.
sel=(s['subset']=='tumor')&(s.dataset=='ALL')&s.signature.isin(['UP','DN'])&(s.adjustment=='unadjusted')&s.P.notna();ii=s.index[sel];pv=s.loc[ii,'P'].to_numpy();order=np.argsort(pv)
if len(ii):
 adj=np.minimum.accumulate((pv[order]*2/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);s.loc[ii[order],'FDR']=adj
save(s,'01_epi/patient_summary.tsv')
print('Verified patient effects',len(pdout),'summary rows',len(s))
