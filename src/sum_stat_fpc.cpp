#include <Rcpp.h>
using namespace Rcpp;

void addToFpc(const NumericMatrix &seg_sites, 
              const NumericVector &positions, 
              NumericMatrix &fpc) {

  bool near;
  bool combinations[2][2] = { false };
  int violations_near = 0,
      violations_far = 0;

  for (int i = 0; i < seg_sites.ncol(); ++i) {
    // Ignore singletons
    if (sum(seg_sites(_,i)) == 1) continue;
    for (int j = i + 1; j < seg_sites.ncol(); ++j) {
      if (sum(seg_sites(_,j)) == 1) continue;
      // Reset combination counter
      combinations[0][0] = false;
      combinations[0][1] = false;
      combinations[1][0] = false;
      combinations[1][1] = false;

      // Are the SNPs near together?
      near = (std::abs(positions[i] - positions[j]) < 0.1); 

      // Count combinations
      for (int k = 0; k < seg_sites.ncol(); ++k) {
        combinations[(int)seg_sites(k,i)][(int)seg_sites(k,j)] = true;
      }

      // If we have all combinations
      if (combinations[0][0] && combinations[0][1] &&
          combinations[1][0] && combinations[1][1]) {
        Rprintf("%i-%i: violates. near: %i\n", i, j, near); 
        if (near) ++violations_near;
        else ++violations_far; 
      }
    }
  }
  ++fpc(violations_near, violations_far);
}

// Exportable wrapper function for unit testing
// [[Rcpp::export]]
NumericMatrix addSegSitesToFpc(const NumericMatrix seg_sites, 
                               const NumericVector positions,
                               NumericMatrix fpc) {

  NumericMatrix fpc_copy = clone(fpc);
  addToFpc(seg_sites, positions, fpc_copy);

  return fpc_copy;
}
