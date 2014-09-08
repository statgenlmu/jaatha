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
  
  try {
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
      size_t digits = 0,  // Number of digits of the length of the tree
             pos = 0,     // Current position on the locus trio
             locus = 0,   // The locus that we are currently in
             len = 0,     // The length of the current tree
             seg_len = 0, 
             locus_end = trio_opts[0]; // The end of the current locus
             
      while (getline(input, line)) {
        // Read in the next tree
        if (line.substr(0, 1) != "[") continue;
        digits = line.find("]")-1;
        len = std::atoi(line.substr(1, digits).c_str());
        
        // If the current tree is valid for a sequence that ends behind the 
        // end of the locus, we need to do a few things:
        while (pos + len >= locus_end) {
          // First print a tree spanning until the end of the current locus 
          // if neccessary.
          seg_len = locus_end - pos;
          //Rprintf("locus %i (end %i) - pos %i - len %i\n", locus, locus_end, 
          //        pos, len);
                  
          if (locus % 2 == 0) {
            output << "[" << seg_len << "]" 
                   << line.substr(digits+2, std::string::npos) << "\n";
          }
          
          // Now the position move towards the end of the current locus.
          len -= seg_len;
          pos += seg_len;
          
          // And look at the next locus.
          if (locus < 4) {
            ++locus;
            locus_end += trio_opts[locus];
          } else {
            // The last tree should end exactly end at the end of the last locus
            if (len != 0) stop("Tree and locus length do not match.");
            // Reset the stats and go one to the next locus trio.
            locus = 0;
            locus_end = trio_opts[0];
            pos = 0;
          }
        }
        
        // If we are here the tree should end within the current locus.
        // Print the rest and move until to the end of the tree.
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
    
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) { 
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  
  return(out_file);
}