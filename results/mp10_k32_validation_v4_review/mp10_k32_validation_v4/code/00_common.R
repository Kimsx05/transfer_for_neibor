.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/refdir/Rlib','/usr/local/lib/R/library'))
suppressPackageStartupMessages({library(data.table);library(qs);library(Matrix);library(digest)})
setDTthreads(2);Sys.setlocale('LC_COLLATE','C');options(stringsAsFactors=FALSE);set.seed(123)
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4'
v3<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3'
old<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1'
out<-function(p)file.path(root,p);wt<-function(x,p)fwrite(x,out(p),sep='\t',na='NA');saveq<-function(x,p){t<-paste0(out(p),'.tmp');qsave(x,t,nthreads=2);stopifnot(file.rename(t,out(p)))}
logmsg<-function(...){cat(format(Sys.time()),...,'\n');flush.console()};memberhash<-function(ids)digest(sort(ids),algo='sha256')
mp<-sprintf('EpiMP%02d',1:14);excluded_sample<-'CNP0000460|CNP0000460_P07T'
config<-list(version='mp10_k32_validation_v4.1',candidate='K32',partition='FULL_Z_r1.00',cluster='C18',n_full=123699L,n_candidate=3282L,n_clusters=27L,seed_new=123L,threads=2L,excluded_sample=excluded_sample,score_controls=list(UCell='C14',NMF='C7'),score_test='two-sided paired Wilcoxon sample means; exact if no ties or zero differences, else normal approximation with continuity correction; Holm over two prespecified tests separately in full and exclude_P07T',DE=yaml::read_yaml(file.path(v3,'CONFIG.yaml'))$DE,fgsea=yaml::read_yaml(file.path(v3,'CONFIG.yaml'))$function_analysis,top_rank='log2FC>0; FDR, PValue, gene_id ascending; no FDR prefilter',patient_analysis='not evaluable unless independently verified mapping exists')
config$fgsea$seed<-123L
if(file.exists(out('CONFIG.yaml')))stopifnot(isTRUE(all.equal(config,yaml::read_yaml(out('CONFIG.yaml')),check.attributes=FALSE)))else yaml::write_yaml(config,out('CONFIG.yaml'),precision=17)
