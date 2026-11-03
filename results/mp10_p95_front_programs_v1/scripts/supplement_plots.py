from plots import *
T=pd.read_csv(O/'tables/summary_effects.tsv',sep='\t');T=T[T.source.eq('tumor_study_specimen_not_region_diagnosis')&T.population.eq('units')]
t=pd.read_csv(O/'tables/ring_summary.tsv.gz',sep='\t');t=t[t.source.eq('tumor_study_specimen_not_region_diagnosis')&t.population_level.eq('units')&t.reference.eq('P95')&t.program.isin(MAIN)]
for pop in ['tissue','Epi']:
 mat=t[t.population.eq(pop)].pivot(index='program',columns='ring',values='mean').reindex(columns=RINGS);norm=mat.sub(mat.mean(axis=1),axis=0).div(mat.std(axis=1).replace(0,np.nan),axis=0);heat(norm,'02_P95_exclusive_rings_'+pop,'Row z-score of mean UCell / '+pop,2)
z=T[T.analysis.isin(['B3_tissue_distance','B3_Epi_distance','B3_tissue_distance_quadratic','B3_Epi_distance_quadratic'])&T.adjustment.eq('adjusted')&T.program.isin(MAIN)];heat(z.pivot(index='program',columns='analysis',values='effect'),'03_P95_distance_effects','Adjusted standardized distance coefficient',.3)
# Region-specific support is separated from the pEMT definition itself.
z=T[T.analysis.eq('A_active_front_support')&T.adjustment.eq('adjusted')&T.program.isin(MAIN)].set_index('program')[['effect']];heat(z,'01_active_front_support','pEMT-high versus other front (adjusted SD)',1)
profile=pd.read_csv(O/'tables/source_R12_profiles.tsv',sep='\t');profile['column']=profile.cluster.astype(str)+' / '+profile.ring
mat=profile.pivot(index='cell_type',columns='column',values='mean_relative_composition');f,a=plt.subplots(figsize=(12,9));im=a.pcolormesh(mat.to_numpy(),cmap='viridis',vmin=0,vmax=.4);a.set_yticks(np.arange(len(mat))+.5,mat.index);a.set_xticks(np.arange(len(mat.columns))+.5,mat.columns,rotation=90);a.invert_yaxis();cb=f.colorbar(im,ax=a);cb.set_label('Original R1/R2 relative composition (display max 0.4)');cb.solids.set_rasterized(False);export(f,'05_original_R12_profiles',[cb.ax])
# Explicit data-set coverage for each major module.
c=pd.read_csv(O/'full_sections.tsv',sep='\t');cov=[]
for ds,g in c.groupby('dataset_id'):
 cov.append(dict(dataset=ds,registered=len(g),Epi=(g.n_query>=20).sum(),P95=(g.n_P95>0).sum(),strict_front=(g.n_front>0).sum(),front_match=(g.n_match_front>=5).sum(),R12=(g.n_R12_query>0).sum()))
v=pd.DataFrame(cov).set_index('dataset');f,a=plt.subplots(figsize=(10,5));im=a.pcolormesh(v.to_numpy(),cmap='Blues');a.set_yticks(np.arange(len(v))+.5,v.index);a.set_xticks(np.arange(len(v.columns))+.5,v.columns,rotation=20,ha='right');a.invert_yaxis()
for i in range(len(v)):
 for j in range(len(v.columns)):a.text(j+.5,i+.5,str(v.iloc[i,j]),ha='center',va='center')
cb=f.colorbar(im,ax=a);cb.set_label('Number of sections');cb.solids.set_rasterized(False);export(f,'00_module_coverage',[cb.ax]);save(v.reset_index(),O/'tables/dataset_module_coverage.tsv')
idx=pd.read_csv(O/'figure_index.tsv',sep='\t');idx=idx[~idx.name.isin([i['name'] for i in index])];save(pd.concat([idx,pd.DataFrame(index)],ignore_index=True),O/'figure_index.tsv')
# Fix vector colorbar export in all section PDFs; PNGs already identical.
os.environ['MP10_PDF_ONLY']='1'
with ProcessPoolExecutor(max_workers=4) as ex:
 for a in ex.map(section_plot,M.section_id):print(a['section_id'],flush=True)
