args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args) && args[1] != "primary") file.path("versions", args[1]) else "."
.libPaths(unique(c(path.expand("~/R/x86_64-pc-linux-gnu-library/4.4"),"/refdir/Rlib",.libPaths())))
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
library(data.table)
setDTthreads(1)
options(mc.cores=1,Ncpus=1)
message("BLAS ", RhpcBLASctl::blas_get_num_procs(), " OMP ", RhpcBLASctl::omp_get_max_threads(), " data.table ", getDTthreads())
skill <- path.expand("~/.agents/skills/jins-routine-plotting")
for(f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Stacked_composition/plot.R")) source(file.path(skill,f))
jin_setup()
cols <- c(ENN="#c95e58",ENS="#bf8257",EIN="#b671a2",EIS="#866aac",NIN="#468b9b",NIS="#56aba2",NNS="#9baf78",NNN="#bcb8b0")
a <- fread(file.path(out,"tables/level_composition_summary.tsv.gz"),cmd=NULL)
a <- a[scope=="overall" & group != "all" & grepl("^region_",measure) & is.finite(effect)]
a[,region:=sub("^region_","",measure)]
a <- a[region %in% names(cols), .(group,region,fraction=effect)]
fwrite(a,file.path(out,"plot_data/architecture_composition.tsv"),sep="\t")
p <- plot_stacked_composition(a,"group","region","fraction",condition_order=intersect(c("infiltrated-like","non-infiltrated-like","mixed","uncertain"),a$group),category_order=names(cols),colors=cols,ylab="Classified tissue proportion")
jin_export(p,file.path(out,"figures/EIS_architecture"),6.3,4.5,legend_plot=p,panel_width=5.1)

