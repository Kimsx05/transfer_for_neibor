
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector score_sparse(IntegerVector ip,IntegerVector ri,NumericVector xx,IntegerVector top,int G,int K,double den){
 int nc=ip.size()-1;NumericVector ans(nc);
 for(int c=0;c<nc;c++){std::vector<double> vals;std::vector<double> tv(top.size(),0);double total=0;
 for(int z=ip[c];z<ip[c+1];z++){double val=xx[z];total+=val;if(val!=0)vals.push_back(val);for(int j=0;j<top.size();j++)if(ri[z]==top[j])tv[j]=val;}
 if(total==0){ans[c]=NA_REAL;continue;}std::sort(vals.begin(),vals.end());double score=0;
 for(double val:tv){double left=std::lower_bound(vals.begin(),vals.end(),val)-vals.begin();double right=std::upper_bound(vals.begin(),vals.end(),val)-vals.begin();double gt=vals.size()-right;double eq=val==0?G-vals.size():right-left;double rank=1+gt+(eq-1)/2;score+=std::max(0.0,K-rank)/den;}ans[c]=score;
 }return ans;}
