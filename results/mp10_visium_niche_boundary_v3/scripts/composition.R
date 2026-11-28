out <- Sys.getenv("V3_MODEL_OUT", ".")
.libPaths(unique(c(path.expand("~/R/x86_64-pc-linux-gnu-library/4.4"),"/refdir/Rlib",.libPaths())))
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
library(data.table)
setDTthreads(1)
options(mc.cores=1,Ncpus=1)
message("BLAS ",RhpcBLASctl::blas_get_num_procs()," OMP ",RhpcBLASctl::omp_get_max_threads()," data.table ",getDTthreads())
skill <- path.expand("~/.agents/skills/jins-routine-plotting")
for(f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Stacked_composition/plot.R")) source(file.path(skill,f))
jin_setup()
cols <- c(ENN="#c95e58",ENS="#bf8257",EIN="#b671a2",EIS="#866aac",NIN="#468b9b",NIS="#56aba2",NNS="#9baf78",NNN="#bcb8b0")
t <- fread(file.path(out,"tables/region_niche_by_section.tsv.gz"))
t <- t[!is.na(unit) & is.finite(fraction) & category!="Unknown"]
s <- t[,.(v=mean(fraction)),by=.(dataset_id,unit,group,type,category)][,.(v=mean(v)),by=.(dataset_id,group,type,category)]
ss <- s[,.(n=uniqueN(group)),by=.(dataset_id,type)][n==3]
a <- merge(s,ss,by=c("dataset_id","type"))[,.(fraction=mean(v),n_sources=uniqueN(dataset_id)),by=.(group,type,category)]
fwrite(a,file.path(out,"tables/composition_plot_data.tsv"),sep="\t")
for(kind in c("region","R12")){
    v <- a[type==kind]
    pal <- if(kind=="region") cols else setNames(c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6"),as.character(0:8))
    v[,category:=as.character(category)]
    p <- plot_stacked_composition(v,"group","category","fraction",condition_order=c("C1","C2","C3"),category_order=names(pal),colors=pal,ylab=if(kind=="region") "Classified tissue proportion" else "Assigned R12 proportion")
    jin_export(p,file.path(out,"figures/groups",paste0("composition_",kind)),5.4,4.7,legend_plot=p,panel_width=4.1)
}

