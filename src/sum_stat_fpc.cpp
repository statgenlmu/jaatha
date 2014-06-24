#include <Rcpp.h>
using namespace Rcpp;

void calcPercentViolation(const NumericMatrix &seg_sites, 
                          const NumericVector &positions,
                          NumericVector &violations, 
                          NumericVector &total_count) {

  bool far;
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
      }
      ++total_count[far];
    }
  }
}



void addToFpc(const NumericMatrix &seg_sites, 
              const NumericVector &positions, 
              const NumericVector &breaks_near,
              const NumericVector &breaks_far,
              NumericMatrix &fpc) {

  NumericVector violations(2); //near - far
  NumericVector total_count(2); //near - far

  calcPercentViolation(seg_sites, positions, violations, total_count);

  // Get percentage of violations inside each class of SNPs
  // and get the locus' class according to between which breaks violations falls.
  NumericVector const *breaks;
  for (int k = 0; k <= 1; ++k) {
    if (k == 0) breaks = &breaks_near;
    else breaks = &breaks_far;

    if (total_count[k] == 0) violations[k] = breaks->size() + 1;
    else {
      violations[k] /= total_count[k];

      for (int i = 0; i < breaks->size(); ++i) {
        if (violations[k] <= (*breaks)[i]) {
          violations[k] = i;
          break;
        }
        if (i == breaks->size() - 1) violations[k] = i + 1;
      }
    }
  }

  ++fpc(violations[0], violations[1]);
}

// Exportable wrapper function for unit testing
// [[Rcpp::export]]
NumericMatrix addSegSitesToFpc(const NumericMatrix seg_sites, 
                               const NumericVector positions,
                               const NumericVector breaks_near,
                               const NumericVector breaks_far,
                               NumericMatrix fpc) {

  NumericMatrix fpc_copy = clone(fpc);
  addToFpc(seg_sites, positions, breaks_near, breaks_far, fpc_copy);

  return fpc_copy;
}

// Exportable wrapper function for unit testing
// [[Rcpp::export]]
NumericVector calcPercentFpcViolation(const NumericMatrix seg_sites, 
                                      const NumericVector positions) {
                             
  NumericVector violations(2);
  NumericVector total_count(2);

  calcPercentViolation(seg_sites, positions, violations, total_count);

  for (int i = 0; i <= 1; ++i) {
    if (total_count(i) == 0) violations(i) = NA_REAL;
    else violations(i) /= total_count(i);
  }

  return violations;
}
