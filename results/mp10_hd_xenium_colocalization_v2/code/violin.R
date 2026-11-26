library(data.table);library(ggplot2);setDTthreads(1)
skill=file.path(path.expand("~"),".agents/skills/jins-routine-plotting")
for (f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Violin_distribution/plot.R"))source(file.path(skill,f))
jin_setup()
for(s in commandArgs(TRUE)){
 d=file.path(if(s %in% c("Cervilla","D1_Johnson","UC_SCC"))"hd" else "xenium",s)
 t=fread(file.path(d,"spatstat_input.tsv.gz"));m0=intersect(names(t),c("MP10","Chol_biosynthesis_8","MP10_covered","H05_transfer"))
 for(m in m0){dt=data.frame(group=ifelse(t[[m]]==1,"high","low"),value=t$density);dt=dt[is.finite(dt$value),]
 p=plot_violin_distribution(dt,"group","value","",order=c("low","high"),colors=c(low="#77a5b0",high="#9d214a"),ylab=if(startsWith(d,"hd"))"Treg component bins/mm²" else "Treg cells/mm²")
 p=p+labs(caption="50 µm · provisional · human review pending");dir.create(file.path("figures",s),showWarnings=FALSE)
 jin_export(p,file.path("figures",s,paste0(m,"_density_violin")),6.5,4.6)
 }
}
