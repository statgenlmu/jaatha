#include <Rcpp.h>
using namespace Rcpp;

size_t idx1, idx2;

void addToPolymClasses(const NumericMatrix seg_sites,
                       const NumericVector sample_size,
                       NumericVector &polym_classes) {
  
  for (int j = 0; j < seg_sites.ncol(); ++j) {
    idx1 = 0;
    idx2 = 0;
    
    for (int i = 0; i < sample_size[0]; ++i) 
      idx1 += seg_sites(i,j); 
    
    for (int i = sample_size[0]; i < seg_sites.nrow(); ++i) 
      idx2 += seg_sites(i,j); 
    
    if (idx1 == 0) {
      if (idx2 == sample_size[1]) ++polym_classes[1];
      else ++polym_classes[0];
    } 
    else if (idx2 == 0) {
      if (idx1 == sample_size[0]) ++polym_classes[1];
      else ++polym_classes[0];
    }
  }
}


NumericVector createPolymVector() {
  return NumericVector::create(
    _["private"] = 0,
    _["fixed"] = 0
  );
}


// For unit testing
// [[Rcpp::export]]
NumericVector classifyPolym(const NumericMatrix seg_sites,
                            const NumericVector sample_size) {
   
  NumericVector classes = createPolymVector();
  addToPolymClasses(seg_sites, sample_size, classes);
  
  return classes;
}