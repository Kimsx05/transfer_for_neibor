
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


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// score_sparse
NumericVector score_sparse(IntegerVector ip, IntegerVector ri, NumericVector xx, IntegerVector top, int G, int K, double den);
RcppExport SEXP sourceCpp_1_score_sparse(SEXP ipSEXP, SEXP riSEXP, SEXP xxSEXP, SEXP topSEXP, SEXP GSEXP, SEXP KSEXP, SEXP denSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ri(riSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type top(topSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type den(denSEXP);
    rcpp_result_gen = Rcpp::wrap(score_sparse(ip, ri, xx, top, G, K, den));
    return rcpp_result_gen;
END_RCPP
}
