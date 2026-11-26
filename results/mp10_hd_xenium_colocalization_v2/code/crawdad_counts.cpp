#include <R.h>
#include <Rinternals.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
// Cache the exact R sample.int(n) order shared by all equal-size tiles at a given offset.
// Geometry and original observation order do not change across label permutations.
extern "C" SEXP tile_counts(SEXP xs, SEXP ys, SEXP marks, SEXP mask, SEXP radius, SEXP reps){
 int n=LENGTH(xs), k=INTEGER(getAttrib(marks,R_DimSymbol))[1], B=asInteger(reps);
 double r=asReal(radius); const double *x=REAL(xs),*y=REAL(ys);
 double x0=*std::min_element(x,x+n), y0=*std::min_element(y,y+n);
 int nx=3+int((*std::max_element(x,x+n)-x0)/r),ny=3+int((*std::max_element(y,y+n)-y0)/r);
 std::vector<std::vector<int>> tiles(nx*ny);
 SEXP out=PROTECT(allocMatrix(REALSXP,B,k));std::fill(REAL(out),REAL(out)+B*k,0.);
 for(int p=0;p<B;++p){
  R_CheckUserInterrupt();for(auto &v:tiles)v.clear();
  for(int i=0;i<n;++i){int a=int((x[i]-x0+p*r/B)/r),b=int((y[i]-y0+p*r/B)/r);tiles[a*ny+b].push_back(i);}
  std::map<int,std::vector<int>> perms;
  for(const auto &v:tiles){
   int m=v.size();if(!m)continue;
   if(std::none_of(v.begin(),v.end(),[&](int i){return LOGICAL(mask)[i];}))continue;
   if(!perms.count(m)){
    SEXP seed=PROTECT(ScalarInteger(p+1));SEXP cs=PROTECT(lang2(install("set.seed"),seed));eval(cs,R_GlobalEnv);
    SEXP sz=PROTECT(ScalarInteger(m));SEXP ca=PROTECT(lang2(install("sample.int"),sz));SEXP pr=PROTECT(eval(ca,R_GlobalEnv));
    perms[m].assign(INTEGER(pr),INTEGER(pr)+m);UNPROTECT(5);
   }
   const auto &perm=perms[m];
   for(int j=0;j<m;++j)if(LOGICAL(mask)[v[j]]){int source=v[perm[j]-1];for(int c=0;c<k;++c)REAL(out)[c*B+p]+=INTEGER(marks)[c*n+source];}
  }
 }
 UNPROTECT(1);return out;
}
