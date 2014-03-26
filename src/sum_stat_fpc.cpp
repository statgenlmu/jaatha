#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix addToFpc(const NumericMatrix seg_sites, NumericMatrix fpc) {
  int c = 0;
  //const NumericVector positions = wrap(seg_sites.attr("dimnames")[1])

  for (int i = 0; i < seg_sites.ncol(); ++i) {
    for (int j = i + 1; j < seg_sites.ncol(); ++j) {
      ++c;
    }
  }

  return fpc;
}
