source('mp10_immune_endothelial_mechanisms_v1/scripts/common.R')
x<-fread(file.path(REF,'collectri_omnipath.tsv'));x[,is_complex:=grepl('COMPLEX',source)];x[,source:=source_genesymbol];x[,target:=target_genesymbol]
# Official decoupleR 2.9.7 get_collectri conventions, frozen fetched OmniPath protein interactions.
x[,mor:=ifelse(is_stimulation==1,1,-1)]
x[is_complex & grepl('JUN|FOS',source),source:='AP1'];x[is_complex & grepl('REL|NFKB',source),source:='NFKB']
x<-unique(x,by=c('source','target'))[,.(source,target,mor)];fwrite(x,file.path(REF,'collectri.tsv'),sep='\t')
e<-new.env();load(file.path(REF,'progeny_model_human_full.rda'),e);p<-as.data.table(e$model_human_full);setorder(p,pathway,p.value,gene);p<-p[,head(.SD,500),by=pathway];setnames(p,c('pathway','gene','weight'),c('source','target','mor'));fwrite(p,file.path(REF,'progeny_top500.tsv'),sep='\t')
e<-new.env();load(file.path(REF,'dorothea_hs.rda'),e);d<-as.data.table(e$dorothea_hs)[confidence %in% c('A','B')];setnames(d,'tf','source');fwrite(d,file.path(REF,'dorothea_AB.tsv'),sep='\t')
lr<-as.data.table(readRDS(file.path(REF,'lr_network_human_21122021.rds')));put(lr,'nichenet_lr_prior.tsv.gz');lt<-readRDS(file.path(REF,'ligand_target_matrix_nsga2r_final.rds'));put(data.table(target=rownames(lt)),'nichenet_target_domain.tsv');put(data.table(ligand=colnames(lt)),'nichenet_ligand_domain.tsv')
print(lr[from %in% c('NECTIN2','PVRL2','PVR','VEGFA') & to %in% c('TIGIT','PVRIG','CD226','FLT1','KDR')]);cat('LT',dim(lt),'\n')
put(CJ(ligand=c('NECTIN2','PVRL2','PVR','VEGFA'),receptor=c('TIGIT','PVRIG','CD226','FLT1','KDR'))[,.(ligand,receptor,in_lr_prior=paste(ligand,receptor) %in% paste(lr$from,lr$to),in_ligand_target_prior=ligand %in% colnames(lt))],'nichenet_key_prior_audit.tsv')
