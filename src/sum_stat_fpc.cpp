#include <Rcpp.h>
#include "seg_sites.h"

using namespace Rcpp;

int findLocus(double position, const NumericVector trio_opts) {
  int ret = 0;
  while (trio_opts[ret] < position) {
    position -= trio_opts[ret];
    ++ret;
  }
  return ret;
}

int far; // 0 = near, 1 = far, 2 = between loci
NumericVector violations; 
NumericVector total_count;

// [[Rcpp::export]]
NumericVector calcPercentFpcViolation(const NumericMatrix seg_sites, 
                                      const NumericVector trio_opts = NumericVector(0)) {
                                        
  NumericVector positions = getPositions(seg_sites);
  
  if (trio_opts.size() == 0) {
    violations = NumericVector(2);
    total_count = NumericVector(2);
  } else {
    violations = NumericVector(3);
    total_count = NumericVector(3);  
  }
  
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
      
      // Overwrite far if they are on different loci
      if (trio_opts.size() > 0 && 
          findLocus(positions[i], trio_opts) != findLocus(positions[j], trio_opts) ) {
          far = 2;
      }

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

  for (int i = 0; i < violations.size(); ++i) {
    if (total_count(i) == 0) violations(i) = NA_REAL;
    else violations(i) /= total_count(i);
  }
  
  return violations;
}

