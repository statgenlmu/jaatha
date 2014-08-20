#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

// Function that reads a file containing output from ms or msms ('in_file')
// and copies the trees in newick format into a separate file ('out_file').
// This new file can be used as input for seq-gen. Only Unix system, this could
// be done with a few grep's, but this function should give us some platform 
// independence.
// [[Rcpp::export]]
std::string parseTrees(std::string in_file, std::string out_file,
                       NumericVector trio_opts) {
    // Open both files
    std::ifstream input(in_file.c_str(), std::ifstream::in);
    std::ofstream output(out_file.c_str(), std::ofstream::out);
    if (!(input.is_open() && output.is_open())) stop("Failed to open file.");
    
    std::string line;
    bool trio = false;
    
    if (trio_opts.size() == 5) {
      Rprintf("%f %f %f %f %f\n", trio_opts[0], trio_opts[1], 
              trio_opts[2], trio_opts[3], trio_opts[4]);
      trio = true;
    }
    
    // Filter trees
    if (!trio) {
      // Fast parser when not using loci-trios.
      while (getline(input, line)) {
        if (line.substr(0, 1) == "[")  output << line << "\n";
      }
    } else {
      // When using loci-trios, we need to remove the inter-locus regions.
      size_t pos = 0, locus = 0, len = 0, locus_end = trio_opts[0];
      while (getline(input, line)) {
        if (line.substr(0, 1) != "[") continue;
        len = std::atoi(line.substr(1, line.find("]")-1).c_str());
        
        if (pos + len == locus_end) {
          // Everything is easy when the tree partition matches the locus' ends
          pos += len;
          ++locus;
          locus_end += trio_opts[locus];
        } else if (pos + len > locus_end) {
          // First move to the end of the current locus / inter-locus region
          
        } else {
          pos += len;
          if (locus % 2 == 1) output << line << "\n";
        }
      }
    }
    
    input.close();
    output.close();
    return(out_file);
}