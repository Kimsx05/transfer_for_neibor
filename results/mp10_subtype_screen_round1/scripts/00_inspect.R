.libPaths(unique(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(qs);library(data.table)})
root <- '/home/data/t070721/codex_workspace/Bladder Metabolism'
out <- file.path(root,'mp10_subtype_screen_round1')
md<-qread(file.path(root,'mp_ucell_state_v1/00_input_audit/cell_metadata.qs'))
cat('EPI METADATA\n');print(names(md));print(head(md[,intersect(names(md),c('cell_id','malignant_like','malignant_like_status','condition','Type','EpiMP10','MP10','cnv_status')),with=FALSE]))
fwrite(md,file.path(out,'tables/epithelial_existing_metadata.tsv.gz'),sep='\t')
cat('MAPPING SUMMARY\n');m<-fread(file.path(root,'subtype/v2/final_annotation_pre_reference/frozen_annotation_barcode_mapping_v2.csv')); print(m[,.N,by=.(major_celltype,annotation_status_v2)]);print(m[,unique(subtype_v2)])
cat('CLINICAL\n');cl<-fread(file.path(root,'cellchat/sample_by_clinical_group.csv'));print(cl[,.N,by=.(dataset,clinical_group,fine_stage,patient_id_status)])
cat('MALIGNANT SUMMARY\n');print(fread('/home/data/t070721/Metabolistic_scRNA/Multi2/3.4 malignant cell process/NMF1/rough_cluster_malignant_like_status_summary.csv'))
cat('SCORE FILES\n'); b<-qread('/home/data/t070721/Metabolistic_scRNA/Multi2/4.3 NMF allgene/output_full_cc_seed5_v2/07_EpiMP_gene_basis.qs');str(b,max.level=1)
suppressPackageStartupMessages(library(CellChat));cat('API\n');print(args(computeCommunProb));print(args(netVisual_circle));print(args(identifyOverExpressedGenes)); print(args(identifyOverExpressedInteractions));print(head(CellChatDB.human$interaction));
