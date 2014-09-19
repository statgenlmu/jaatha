#include <Rcpp.h>
using namespace Rcpp;

bool far;

// [[Rcpp::export]]
NumericVector calcPercentFpcViolation(const NumericMatrix seg_sites, 
                                      const NumericVector positions) {
                                        
                             
  NumericVector violations(2);
  NumericVector total_count(2);
  bool combinations[2][2] = { false };

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
      far = (std::abs(positions[i] - positions[j]) > 0.1); 

      // Count combinations
      for (int k = 0; k < seg_sites.nrow(); ++k) {
        combinations[(int)seg_sites(k,i)][(int)seg_sites(k,j)] = true;
      }

      // If we have all combinations
      if (combinations[0][0] && combinations[0][1] &&
          combinations[1][0] && combinations[1][1]) {
        ++violations[far];
        //Rprintf("%i %i: far:%i vio:1\n", i, j, far);
      } else {
        //Rprintf("%i %i: far:%i vio:0\n", i, j, far);
      }
      ++total_count[far];
    }
  }

  for (int i = 0; i <= 1; ++i) {
    if (total_count(i) == 0) violations(i) = NA_REAL;
    else violations(i) /= total_count(i);
  }

  return violations;
}
