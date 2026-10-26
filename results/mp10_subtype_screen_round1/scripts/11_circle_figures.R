.libPaths(c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(qs);library(CellChat)})
root<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1';tb<-function(n)file.path(root,'tables',n)
axes<-fread(tb('selected_axes.tsv'));ed<-fread(tb('CellChat_detected_edges.tsv.gz'));pal<-fread(tb('cell_colors.tsv'));colors<-setNames(pal$color,pal$cell_type);db<-fread(tb('CellChat_DB_interactions.tsv'));files<-list.files(file.path(root,'objects/networks'),full.names=TRUE,pattern='qs$');nets<-setNames(lapply(files,qread),sub('.qs$','',basename(files)))
specs<-list(Normal_Tumor=c('Normal','Tumor'),T1_T2=c('T1','T2'),MP10_low_high=c('low','high'))
safe<-function(z)gsub('[^A-Za-z0-9_-]','_',z);erows<-list();nrows<-list();ii<-0L;jj<-0L
edge_value<-function(net,src,tgt,lr){
 if(is.null(net))return(list(pr=NA_real_,pv=NA_real_,value=NA_real_,status='missing_condition'))
 dn<-dimnames(net$net$prob)
 if(!src%in%dn[[1]]||!tgt%in%dn[[2]])return(list(pr=NA_real_,pv=NA_real_,value=NA_real_,status='source_or_target_below_30'))
 if(!lr%in%dn[[3]])return(list(pr=NA_real_,pv=NA_real_,value=0,status='not_detected_LR_screen'))
 pr<-net$net$prob[src,tgt,lr];pv<-net$net$pval[src,tgt,lr];list(pr=pr,pv=pv,value=if(pv<.05)pr else 0,status=if(pv<.05&&pr>0)'detected' else 'not_detected_probability')
}
for(ai in seq_len(nrow(axes))){a<-axes[ai];
 # Fixed display nodes are chosen once for this LR and sender across all scopes.
 extra<-ed[source==a$source & LR_ID==a$LR_ID & !target%in%c(a$source,a$target),.(w=max(probability)),by=target][order(-w)]
 nodes<-unique(c(a$source,a$target,head(extra$target,4)));nn<-length(nodes);layout<-igraph::norm_coords(igraph::layout_in_circle(igraph::make_empty_graph(nn,directed=TRUE)),xmin=-1,xmax=1,ymin=-1,ymax=1);rownames(layout)<-nodes
 for(ds in c('ALL',sort(unique(ed[scope!='ALL',scope]))))for(co in names(specs)){
  gs<-specs[[co]];ns<-lapply(gs,function(g)nets[[paste(ds,co,g,sep='__')]])
  if(all(vapply(ns,is.null,logical(1))))next
  mats<-list();states<-list();localrows<-list()
  for(k in 1:2){g<-gs[k];net<-ns[[k]];mat<-matrix(0,nn,nn,dimnames=list(nodes,nodes));rr<-list()
   for(tgt in setdiff(nodes,a$source)){
    z<-edge_value(net,a$source,tgt,a$LR_ID);if(is.finite(z$value))mat[a$source,tgt]<-z$value
    rr[[tgt]]<-data.table(candidate=a$candidate_subtype,source=a$source,target=tgt,LR_ID=a$LR_ID,LR_display=a$LR_display,pathway=a$pathway,scope=ds,comparison=co,condition=g,probability=z$pr,internal_P=z$pv,plotted_strength=z$value,status=z$status)
   }
   localrows[[k]]<-rbindlist(rr);mats[[k]]<-mat
   states[[k]]<-if(all(is.na(localrows[[k]]$plotted_strength)))'Not evaluable' else if(sum(mat)==0)'Not detected under current criteria' else ''
   for(v in nodes){jj<-jj+1L;nc<-if(is.null(net))0 else net$counts[communication_type==v,N];nrows[[jj]]<-data.table(candidate=a$candidate_subtype,LR_ID=a$LR_ID,scope=ds,comparison=co,condition=g,node=v,order=match(v,nodes),x=layout[v,1],y=layout[v,2],color=colors[v],node_weight=1,n_cells=if(length(nc))nc else 0,evaluable=length(nc)>0&&nc>=30)}
  }
  mx<-max(vapply(mats,max,numeric(1)));sc<-max(mx,1e-12)
  for(k in 1:2){ii<-ii+1L;erows[[ii]]<-localrows[[k]][,`:=`(edge_weight_max=mx,edge_width_max=8,node_set=paste(nodes,collapse=';'))]}
  folder<-if(ds=='ALL')file.path(root,'figures/overall') else file.path(root,'figures/by_dataset',ds);dir.create(folder,recursive=TRUE,showWarnings=FALSE);stem<-file.path(folder,paste0('circle_',co,'_',safe(a$candidate_subtype),'_',a$LR_ID))
  draw<-function(){
   par(mfrow=c(1,2),mar=c(3,2,4,2),oma=c(1,0,3,0),xpd=NA)
   for(k in 1:2){
    netnow<-ns[[k]];present<-if(is.null(netnow))character() else netnow$counts[N>=30,communication_type]
    displaymat<-mats[[k]];displaylabels<-ifelse(nodes%in%present,nodes,paste0(nodes,'\n[NE]'));dimnames(displaymat)<-list(displaylabels,displaylabels)
    CellChat::netVisual_circle(displaymat,color.use=unname(colors[nodes]),title.name=gs[k],layout=igraph::in_circle(),vertex.weight=rep(1,nn),vertex.weight.max=1,vertex.size.max=12,vertex.label.cex=.78,vertex.label.color="#00000000",edge.weight.max=sc,edge.width.max=8,weight.scale=TRUE,remove.isolate=FALSE,arrow.size=.45,arrow.width=1.3,margin=.4)
    # Position labels outside nodes explicitly; native CellChat labels can overlap long names.
    labels_draw<-vapply(displaylabels,function(z)paste(strwrap(gsub('Malignant_like_epithelial','Malignant_like_\nepithelial',z,fixed=TRUE),width=24),collapse='\n'),character(1))
    for(vi in seq_len(nn)){xx<-layout[vi,1];yy<-layout[vi,2];text(xx+sign(xx)*.14,yy+sign(yy)*.14,labels_draw[vi],adj=c(if(xx>.1)0 else if(xx<(-.1))1 else .5,.5),cex=.74,xpd=NA)}
    if(nzchar(states[[k]]))mtext(states[[k]],side=1,line=.4,cex=.85)
   }
   mtext(a$LR_display,side=3,outer=TRUE,line=1,cex=1.25);mtext(paste(ds,co,'| sender:',a$source),side=3,outer=TRUE,line=-.5,cex=.8)
   mtext(paste('Inferred LR direction; equal nodes; shared edge maximum =',format(mx,digits=3)),side=1,outer=TRUE,line=-.3,cex=.78)
  }
  pdf(paste0(stem,'.pdf'),width=14,height=6.5,bg='transparent',useDingbats=FALSE);draw();dev.off()
  png(paste0(stem,'.png'),width=14,height=6.5,units='in',res=300,bg='transparent');draw();dev.off()
  if(ds=='ALL'&&co=='MP10_low_high'){png(file.path(root,'logs',paste0('preview_white_',safe(a$candidate_subtype),'.png')),width=1960,height=910,res=140,bg='white');draw();dev.off()}
 }
}
fwrite(rbindlist(erows),tb('circle_edges.tsv'),sep='\t',na='NA');fwrite(rbindlist(nrows),tb('circle_nodes.tsv'),sep='\t',na='NA');cat('Circle figures complete\n')
