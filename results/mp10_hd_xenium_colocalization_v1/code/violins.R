source(file.path(Sys.getenv("OUT"),"code/io.R"))
library(ggplot2)
skill<-file.path(path.expand("~"),".agents/skills/jins-routine-plotting")
for(f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Violin_distribution/plot.R"))source(file.path(skill,f))
jin_setup()
for(s in c("A1R2","A1R4","A1R5","A1R6","A2R2")){
 d<-file.path(out,"xenium",s);scores<-fread(file.path(d,"ucell.tsv.gz"));v<-melt(scores,id.vars="id",measure.vars=c("Chol_biosynthesis_8","MP10_covered"),variable.name="program",value.name="score")
 p<-plot_violin_distribution(v,"program","score","UCell",order=c("Chol_biosynthesis_8","MP10_covered"),colors=c(Chol_biosynthesis_8="#9c3f67",MP10_covered="#487da3"),ylab="Score")
 jin_export(p,file.path(out,"figures",s,"measured_state_distribution"),width=6.5,height=4.6)
}

