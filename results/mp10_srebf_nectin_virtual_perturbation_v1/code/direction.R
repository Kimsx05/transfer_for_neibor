# Exact Pearson Gram matrix for genes, using BLAS to avoid pairwise scalar loops.
gram_direction=function(ct,genes,ko,depth=FALSE) {
 S=scTenifoldKnk:::.logCPM(ct,genes)
 L=if(depth)scTenifoldKnk:::.regressLibSize(S) else S$L-S$mu
 L=L-rowMeans(L);ss=sqrt(rowSums(L^2));ss[ss==0]=Inf
 R=tcrossprod(L/ss);R[!is.finite(R)]=0
 H=scTenifoldKnk::heatKernel(R,t=5)
 ks=as.list(ko);names(ks)=ko
 scTenifoldKnk:::.kernelEffect(H,S,ks)
}

