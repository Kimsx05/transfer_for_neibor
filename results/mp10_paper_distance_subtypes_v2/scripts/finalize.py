from plots import *
import zipfile,shutil
def main():
 sc=pd.read_csv(O/'tables/sc_sample_programs.tsv',sep='\t');heat(sc.groupby('type')[PROG].median().T,'Single-cell programs / median of sample medians','single_cell',.3)
 re=pd.read_csv(O/'tables/reactome_enrichment.tsv',sep='\t');v=re[re.dataset.eq('dataset_equal')].sort_values('padj').head(25).set_index('pathway')
 heat(v[['NES']],'Reactome / top 25 adjusted-MP10 rank enrichments (full tests retained)','reactome',3)
 a=pd.read_csv(O/'tables/summary.tsv',sep='\t');a=a[a.source.str.startswith('tumor')&a.population.eq('units')]
 t=pd.read_csv(O/'tables/targets.tsv.gz',sep='\t');rows=[]
 for name in list(TARGETS)[:7]:
  e=a[(a.kind=='distance')&a.axis.eq(name+'_P95')];mp=e[e.program==RAW].iloc[0];fun=e[e.program!=RAW].sort_values('q').head(3);cov=t[(t.target==name)&t['quantile'].eq(95)]
  rows.append(dict(subtype=name,sections_defined=cov.n_target.gt(0).sum(),median_target_abundance=cov.abundance_median.median(),low_abundance_sections=cov.low_abundance.sum(),adjusted_MP10_r=mp.effect,CI_low=mp.low,CI_high=mp.high,q=mp.q,positive_datasets=mp.positive_datasets,n_datasets=mp.n_datasets,programs=' | '.join(f'{r.program} r={r.effect:.3f} q={r.q:.3g}' for r in fun.itertuples()),spatial_context='section / v3.4 / R12 / local tables; local interface coverage limited',hypothesis='subtype-specific program and location association for follow-up',missing_evidence='independent replication; cell-specific paired state; perturbation; no robust main interaction'))
 write(rows,O/'tables/candidates.tsv')
 for f in ['gene_sets.tsv.gz','gene_set_coverage.tsv','source_R12_profiles.tsv']:
  shutil.copyfile(F/'tables'/f,O/'tables'/f)
 shutil.copyfile(F/'tables/associations.tsv',O/'tables/historical_P80_reference.tsv')
 figs=[]
 for p in sorted((O/'figures').glob('*.png')):
  figs.append(dict(figure=p.stem,pdf=str(p.with_suffix('.pdf').relative_to(O)),png=str(p.relative_to(O))))
 write(figs,O/'figure_index.tsv')
 for ds,g in M.groupby('dataset_id'):
  files=[O/'figures'/f'{s}_{kind}.pdf' for s in g.section_id for kind in ['targets','interface']]
  subprocess.run(['Rscript','-e','a=commandArgs(TRUE);qpdf::pdf_combine(a[-1],a[1])',str(O/f'{ds}_atlas.pdf'),*map(str,files)],check=True)
 names=['sensitivity_comparison','interface_composition','target_support','MP10_effects','program_distance','MP10_programs','interactions','interface_programs','interface_niche','niche_interactions','single_cell','reactome']+['curves_'+str(i) for i in range(1,8)]
 subprocess.run(['Rscript','-e','a=commandArgs(TRUE);qpdf::pdf_combine(a[-1],a[1])',str(O/'overview.pdf'),*[str(O/'figures'/f'{n}.pdf') for n in names]],check=True)
 with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED) as z:
  files=[O/'overview.pdf',O/'METHODS.md',O/'RESULTS.md',O/'figure_index.tsv']
  files += [p for p in (O/'tables').glob('*') if not any(p.name.startswith(s+'_') for s in M.section_id)]
  for p in files:z.write(p,p.relative_to(O))
if __name__=='__main__':main()
