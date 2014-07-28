#include <Rcpp.h>
using namespace Rcpp;

// Fast function to add a locus' seg.sites to a JSFS
void addToJsfs(const NumericMatrix &seg_sites,
               const NumericVector &sample_size,
               NumericMatrix &jsfs) {

  int idx1;
  int idx2;
  
  for (int j = 0; j < seg_sites.ncol(); ++j) {
    idx1 = 0;
    idx2 = 0;

    for (int i = 0; i < sample_size[0]; ++i) idx1 += seg_sites(i,j); 
    for (int i = sample_size[0]; i < seg_sites.nrow(); ++i) idx2 += seg_sites(i,j); 

    ++jsfs(idx1, idx2);
  }
}

// Exportable wrapper function for unit testing
// [[Rcpp::export]]
NumericMatrix addSegSitesToJsfs(const NumericMatrix seg_sites, 
                                const NumericVector sample_size,
                                NumericMatrix jsfs) {

  NumericMatrix jsfs_copy = clone(jsfs);
  addToJsfs(seg_sites, sample_size, jsfs_copy);

  return jsfs_copy;
}
