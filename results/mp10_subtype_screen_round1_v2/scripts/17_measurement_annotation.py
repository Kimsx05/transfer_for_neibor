from common import *
c=read('signaling_gene_dataset_coverage.tsv');available={ds:set(z.loc[z.observed_in_dataset,'gene']) for ds,z in c.groupby('dataset')};a=read('candidate_sample_expression.tsv.gz');a['genes_observed_somewhere_in_dataset']=[set((x.ligand_genes+';'+x.receptor_genes).split(';')).issubset(available[x.dataset]) for x in a.itertuples()];a['gene_measurement_status']=np.where(a.genes_observed_somewhere_in_dataset,'dataset_nonzero_evidence; node-specific zeros retained','feature_present_but_all_zero_in_dataset; measurement_unresolved')
# Never interpret a cohort-wide structural all-zero feature as demonstrated measured absence.
a.loc[~a.genes_observed_somewhere_in_dataset,['epi_delta','aux_product_delta']]=np.nan
save(a,'candidate_sample_expression.tsv.gz')
b=read('B_marker_sample_support.tsv');bc=b.groupby(['dataset','gene']).agg(n_B_cells=('n_cells','sum'),n_specimen_subtype_groups=('n_cells','size'),observed_nonzero_in_B=('mean',lambda x:(x>0).any())).reset_index();bc['feature_present_in_frozen_RNA']=bc.gene.isin(read('B_marker_gene_coverage.tsv').query('measured').gene);bc['original_dataset_feature_manifest']='not_available; all-zero in B does not demonstrate unmeasured gene';save(bc,'B_marker_dataset_gene_coverage.tsv')
print('Measurement unresolved rows',(~a.genes_observed_somewhere_in_dataset).sum())
