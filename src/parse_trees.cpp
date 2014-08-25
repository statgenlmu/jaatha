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
    
    if (trio_opts.size() == 5) trio = true;
    
    // Filter trees
    if (!trio) {
      // Fast parser when not using loci-trios.
      while (getline(input, line)) {
        if (line.substr(0, 1) == "[")  output << line << "\n";
      }
    } else {
      // When using loci-trios, we need to remove the inter-locus regions.
      size_t digits = 0, pos = 0, locus = 0, len = 0, seg_len = 0, 
             locus_end = trio_opts[0];
             
      while (getline(input, line)) {
        if (line.substr(0, 1) != "[") continue;
        digits = line.find("]")-1;
        len = std::atoi(line.substr(1, digits).c_str());
        
        while (pos + len >= locus_end) {
          // Print tree until end of locus if neccessary
          seg_len = locus_end - pos;
          //Rprintf("locus %i (end %i) - pos %i - len %i\n", locus, locus_end, 
          //        pos, len);
                  
          if (locus % 2 == 0) {
            output << "[" << seg_len << "]" 
                   << line.substr(digits+2, std::string::npos) << "\n";
          }
          
          // Move to the end of the locus
          len -= seg_len;
          pos += seg_len;
          ++locus;
          locus_end += trio_opts[locus];
          
          // Reset counter if we reached end of last locus.
          if (locus == 5 && len == 0) {
            locus = 0;
            locus_end = trio_opts[0];
            pos = 0;
            len = 0;
          }
        }
        
        // If length is left, this now falls completely inside the current 
        // locus.
        if (len > 0) {
           pos += len;
           if (locus % 2 == 0) {
            output << "[" << len << "]" 
                   << line.substr(digits+2, std::string::npos) << "\n";
           }
        }
      }
      
      if (pos != 0) stop("Error parsing trees");
    }
    
    input.close();
    output.close();
    return(out_file);
}