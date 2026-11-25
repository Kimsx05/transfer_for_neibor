from common import *
T=Path('tables');O=Path('plot_data');read=lambda f:rd(T/f)
cats=read('rank_category.tsv');v=cats[cats.analysis.isin(['v1','full','clean']) & cats.ko.isin(K) & (cats.branch=='heat') & (cats.top_n==50)].copy();v['category']=v.category.replace({'mitochondrial_genome':'MT genome','mitochondrial_ribosome':'Mitochondrial ribosome','cytosolic_ribosome':'Cytosolic ribosome','other':'Other','unclassified':'Unclassified'});save(v,O/'F01.tsv')
a=read('control_specificity.tsv');g=a[a['rank']<=8].gene.unique();save(a[a.gene.isin(g)],O/'F02.tsv')
save(a[a.gene.isin(['NECTIN1','NECTIN2','NECTIN3','NECTIN4','PVR','FASN','SCD','MARCHF6','RABGEF1','HMGCR','SQLE','RPS6KB1','LPIN1','EIF4EBP1'])],O/'F03.tsv')
p=read('mtor_module_summary.tsv');save(p[(p.kind=='program')&(p['mask']=='Full')],O/'F04.tsv')
r=read('robustness.tsv');save(r[(r.gene=='NECTIN2')&r.refit_type.isin(['LODO','patient_subsample','state','seed','clean'])],O/'F05.tsv')
save(a[a.gene.isin(['NECTIN2','FASN','SCD','MARCHF6','RABGEF1','HMGCR','PRDX1'])],O/'F06.tsv')
c=read('cross_method_summary.tsv');save(c[(c.analysis=='main')&c.ko.isin(K[:2]) & c.gene.isin(['NECTIN1','NECTIN2','NECTIN3','NECTIN4','PVR','ACACA','FASN','HMGCR','SCD','SQLE','LDLR','INSIG1','MARCHF6'])],O/'F07.tsv')
pair=read('patient_paired_summary.tsv');save(pair[pair.gene=='NECTIN2'],O/'F08.tsv')

