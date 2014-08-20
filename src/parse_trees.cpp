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
    
    while (getline(input, line)) {
      if (line.substr(0, 1) == "[") output << line << "\n";
    }
    
    input.close();
    output.close();
    return(out_file);
}