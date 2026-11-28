from common3 import *
z=read(O/'cache/immune_features.tsv');md=[]
for (grid,line),g in z.merge(M.tumor.rename('tumor'),left_on='section_id',right_index=True).groupby(['grid','axis']):
    v=g[g.tumor]
    for k in list('ABCDEFGP'):
        md.append(dict(grid=grid,axis=line,feature=k,tumor_sections=len(v),finite=v[k].notna().sum(),true_zero=v[k].eq(0).sum(),zero_variance=v[k].std(ddof=0)==0,missing_reason='none' if v[k].notna().all() else 'insufficient eligible windows, no actual E-neighbor pairs/variance or no evaluable outer domain; inspect per-section support'))
save(md,O/'audit/candidate_coverage.tsv')
w=pd.read_parquet(O/'cache/immune_windows.parquet').query("grid=='500' and eligible");w2=read(V/'classification/immune_windows.tsv.gz').query("level=='window500'")[['section_id','local_id','per_n']].rename(columns={'local_id':'window','per_n':'v2_same_grid_per_n'})
w=w.merge(w2,on=['section_id','window'],validate='1:1');tmp=[]
for s,a in w.groupby('section_id'):
    def wm(k):
        v=a.sort_values(k);return v[k].to_numpy()[np.searchsorted(v.E_n.cumsum(),v.E_n.sum()/2)]
    high=a.frac_T>=wm('frac_T');dense=a.density_T>=wm('density_T')
    tmp.append(pd.DataFrame(dict(section_id=s,window=a.window,E_n=a.E_n,frac_T=a.frac_T,density_T=a.density_T,above_section_fraction_median=high,above_section_density_median=dense,discordant=high!=dense,assigned_crossgrid_outer_n=a.outer_n,v2_own_grid_outer_n=a.v2_same_grid_per_n)))
save(pd.concat(tmp),O/'audit/window_density_fraction_and_outer.tsv.gz')
ev=read(V/'audit/section_clinical.tsv').merge(read(O/'audit/literature_key.tsv'),on='dataset_id').merge(M.image_path.apply(lambda f:Path(f).is_file()).rename('HE_exists'),left_on='section_id',right_index=True)
ev['original_article']=ev.doi.map(lambda s:'https://doi.org/'+s if pd.notna(s) and s!='unknown' else None);save(ev,O/'audit/section_literature.tsv')

