suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrastr);library(jsonlite)})
a<-commandArgs(TRUE);root<-normalizePath(a[1]);s<-a[2];out<-file.path(root,s);fig<-file.path(out,'figures');dir.create(fig,showWarnings=FALSE)
skill<-'/home/data/t070721/.agents/skills/jins-routine-plotting';for(f in c('theme','palettes','export'))source(file.path(skill,'common',paste0(f,'.R')));jin_setup()
coarse<-c(Epithelial='#f5b264',Mesenchymal='#8b78b8',Endothelial='#2878b5',Immune='#57c1b5',Unknown='#a8a8a8')
export<-function(p,name,w=9,h=6,pw=w){if(file.exists(file.path(fig,paste0(name,'.png'))))return(invisible());jin_export(p,file.path(fig,name),width=w,height=h,panel_width=pw)}
spfinish<-function(p)p+coord_fixed()+jin_theme()+labs(x=if(grepl('MPBC',s))'Relative array x (um)' else 'Original x (um)',y=if(grepl('MPBC',s))'- Relative array y (um)' else '- Original y (um)')
roi_fixed<-function(p,b){p+coord_fixed(xlim=c(b$x,b$X),ylim=c(-b$Y,-b$y),expand=FALSE)+annotate('segment',x=b$x+12,xend=b$x+76,y=-b$Y+12,yend=-b$Y+12,linewidth=1.4,color='white')+annotate('segment',x=b$x+12,xend=b$x+76,y=-b$Y+12,yend=-b$Y+12,linewidth=.4,color='black')+annotate('label',x=b$x+44,y=-b$Y+28,label='64 um',fill='white',linewidth=0,size=3)}
