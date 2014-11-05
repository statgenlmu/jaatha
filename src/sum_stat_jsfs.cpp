#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the JSFS from a list of segregating sites statistics
// [[Rcpp::export]]
NumericMatrix calcJsfs(const List seg_sites, const NumericVector sample_size) {
  
  NumericMatrix jsfs(sample_size[0]+1, sample_size[1]+1);
  size_t idx1, idx2;
  
  for (int locus = 0; locus < seg_sites.size(); ++locus) {
    NumericMatrix ss = as<NumericMatrix>(seg_sites[locus]);
    for (int j = 0; j < ss.ncol(); ++j) {
      idx1 = 0;
      idx2 = 0;
    
      for (int i = 0; i < sample_size[0]; ++i) idx1 += ss(i,j); 
      for (int i = sample_size[0]; i < ss.nrow(); ++i) idx2 += ss(i,j); 
    
    ++jsfs(idx1, idx2);
    }
  }
  
  return jsfs;
}