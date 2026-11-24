source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
# Read compact evidence and restrict the human v2 ligand-target prior to M0073.
lt=readRDS(file.path(R,"ligand_target_matrix_nsga2r_final.rds"));d=tab(B,"de/M0073.tsv.gz");bg=intersect(d$gene,rownames(lt))
ids=tab(V,"mTORC1_M0073_membership.tsv")$sample_uid
det=tab(B,"prior_molecules_per_profile.tsv.gz")[sample_uid%in%ids & n_cells>=50 & !profile_type%in%c("H50","H10","H05","L50"),.(n_profiles=.N,fraction=median(positive_fraction)),by=.(gene,profile_type)]
lr=unique(tab(B,"nichenet_lr_prior.tsv.gz")[,.(ligand=from,receptor=to)])
rec=det[profile_type=="Capillary EC" & fraction>=.10,.(receptor=gene,receptor_detection=fraction)]
lr=merge(lr,rec,by="receptor");lr=lr[ligand%in%colnames(lt)]
sources=det[gene%in%lr$ligand&fraction>=.1];put(sources,"communication_sender_detection.tsv")
edges=merge(lr,sources[profile_type=="Epi_all",.(ligand=gene,epi_detection=fraction)],by="ligand",all.x=TRUE)
le=gettab("legacy_full_m20.tsv")[mask=="Full" & pathway%in%paste0("HALLMARK_",c("MYC_TARGETS_V1","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","MTORC1_SIGNALING"))]
pr=lt[bg,unique(lr$ligand),drop=FALSE];rm(lt);gc()
activity=list()
for(k in seq_len(nrow(le))){
 targets=intersect(strsplit(le$leadingEdge[k],";")[[1]],bg);y=bg%in%targets
 met=rbindlist(lapply(seq_len(ncol(pr)),function(j){x=pr[,j];if(diff(range(x))==0)return(data.table(ligand=colnames(pr)[j],aupr_corrected=NA_real_,pearson=NA_real_))
 pred=ROCR::prediction(x,y);ro=ROCR::performance(pred,"prec","rec");precision=ro@y.values[[1]];precision[is.na(precision)]=1
 data.table(ligand=colnames(pr)[j],aupr_corrected=caTools::trapz(ro@x.values[[1]],precision)-mean(y),pearson=cor(x,as.integer(y)))}))
 met[,program:=le$pathway[k]];met[,n_targets:=length(targets)];met[,n_background:=length(bg)];met[,epi_qualified:=ligand%in%edges[!is.na(epi_detection),ligand]];setorder(met,-aupr_corrected,ligand,na.last=TRUE)
 met[is.finite(aupr_corrected),agnostic_rank:=seq_len(.N)];met[epi_qualified & is.finite(aupr_corrected),epi_rank:=seq_len(.N)]
 activity[[k]]=met
}
act=rbindlist(activity);put(act,"leading_edge_ligand_ranking.tsv.gz")
links=merge(act[epi_qualified==TRUE],edges,by="ligand",allow.cartesian=TRUE);put(links,"leading_edge_LR.tsv.gz")
cc=rbindlist(lapply(Sys.glob("mp10_cellchat_high_definition_sensitivity_v1/tables/lr/C10__ALL__*__raw_sparse.tsv.gz"),function(f)fread(f)[target=="Capillary EC"&supported==TRUE]))
put(cc,"CellChat_capillary_both_directions.tsv.gz")
rev=tab(B,"nichenet_ligand_activities_annotated.tsv.gz")[model_id%in%c("M0001","M0016","M0031") & mask=="M20-masked" & detection_threshold==.1 &target_direction=="positive"&grepl("Capillary EC",detected_sender_types)]
put(rev,"reverse_capillary_sources.tsv")
put(fread("mp10_spatial_communication_v1/epi_state_type_contrasts.tsv")[lineage=="Endothelial"],"spatial_bidirectional.tsv")
put(tab(V,"mTORC1_spatial_source_evidence.tsv"),"spatial_mTOR_source.tsv")

