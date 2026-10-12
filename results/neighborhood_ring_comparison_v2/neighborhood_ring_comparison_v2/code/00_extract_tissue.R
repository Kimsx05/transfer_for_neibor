args<-commandArgs(TRUE)
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(Seurat))
x<-qread(args[1],nthreads=2)
y<-list(); j<-0L
for(i in seq_along(x)) {
 s<-names(x)[i]
 if(!s %in% names(x[[i]]@images)) next
 co<-x[[i]]@images[[s]]@coordinates
 j<-j+1L
 y[[j]]<-data.frame(section=s,barcode=rownames(co),array_row=co$row,array_col=co$col,tissue=if('tissue'%in%colnames(co)) co$tissue else NA,imagerow=co$imagerow,imagecol=co$imagecol)
}
write.csv(do.call(rbind,y),gzfile(args[2]),row.names=FALSE)
