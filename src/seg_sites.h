#ifndef jaatha_src_seg_sites
#define jaatha_src_seg_sites

namespace Rcpp {
  
inline NumericVector getPositions(NumericMatrix seg_sites) {
  if (!seg_sites.hasAttribute("positions")) stop("SegSites without positions");
  return seg_sites.attr("positions");
}
  
}

#endif