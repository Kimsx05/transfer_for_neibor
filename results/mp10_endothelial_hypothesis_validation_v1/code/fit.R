fitone=function(x,y,dataset,cov=NULL){
 d=data.frame(x=x,y=y,ds=factor(dataset));if(!is.null(cov))d$cov=cov
 d=d[complete.cases(d),,drop=FALSE];d$ds=droplevels(d$ds);n=nrow(d)
 ret=data.table(n,df=NA_integer_,beta=NA_real_,lo=NA_real_,hi=NA_real_,P=NA_real_,partial_r=NA_real_,partial_rho=NA_real_,status="NE coverage or rank")
 if(n<4)return(ret)
 adj=c(if(nlevels(d$ds)>1)"ds",if(!is.null(cov))"cov")
 fit=lm(reformulate(c(adj,"x"),"y"),d);ret[,df:=fit$df.residual]
 if(fit$df.residual<3||fit$rank!=length(coef(fit))||sd(d$y)==0)return(ret)
 z=coef(summary(fit))["x",];ci=confint(fit,"x")
 res=function(a){if(!length(adj))return(a-mean(a));resid(lm(reformulate(adj,"a"),data=transform(d,a=a)))}
 ret[, `:=`(beta=z[1],lo=ci[1],hi=ci[2],P=z[4],partial_r=cor(res(d$x),res(d$y)),partial_rho=cor(res(rank(d$x)),res(rank(d$y))),status="completed")];ret
}
