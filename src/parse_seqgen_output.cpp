#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;


NumericMatrix parseSeqgenSegSites(std::ifstream &output,
                                  size_t locus_length,
                                  const size_t individuals,
                                  NumericVector &position,
                                  const NumericVector trio_opts) {
  
  std::string tmp;
  size_t seq_nr;
  
  // First read the complete locus and save it in `sequence`.
  std::vector<std::vector<char> > sequence(individuals);
  for (size_t i=0; i<individuals; ++i) {
    // Read sequence number
    if (!output.good()) 
    stop("Unexpeced end of seq-gen file.");
    output >> tmp;
    // ms from phyclust adds an "s" to the seqName. Remove it if there.
    if (tmp.compare(0, 1, "s") == 0) tmp.erase(0,1);
    seq_nr = atoi(tmp.c_str());
    
    // Read sequence
    output >> tmp;
    const char* cstr=tmp.c_str();
    sequence.at(seq_nr-1).assign(cstr, cstr+tmp.length());
  }
  
  // Determine which positions are SNPs
  std::vector<double> positions;
  size_t derived_count;
  
  for (size_t j=0; j<locus_length; ++j) {
    derived_count = 0;
    for (size_t i=0; i<individuals-1; ++i) {
      derived_count += (sequence.at(i).at(j) != sequence.at(individuals-1).at(j));
    }
    if (derived_count > 0 && derived_count < (individuals - 1)) {
      positions.push_back(j);
    }
  }
  
  NumericMatrix seg_sites(individuals-1, positions.size());
  
  if (positions.size() > 0) {
    for (size_t i=0; i<individuals-1; ++i) {
      derived_count = 0;
      for (std::vector<double>::iterator it = positions.begin(); it != positions.end(); ++it) {
        seg_sites(i, derived_count) = (sequence[i][*it] != sequence[individuals-1][*it]);
        ++derived_count;
      }
    }
    
    if (trio_opts.size() != 5) {
      for (std::vector<double>::iterator it = positions.begin(); 
      it != positions.end(); ++it) {
        *it /= (locus_length - 1);
      }
    } else {
      locus_length = sum(trio_opts);
      size_t offset = 0, locus = 0, locus_end = trio_opts[0];
      for (std::vector<double>::iterator it = positions.begin(); 
      it != positions.end(); ++it) {
        while(*it >= locus_end) {
          ++locus;
          offset += trio_opts[locus];   // Increase offset that locus 1 + 3
          ++locus;
          locus_end += trio_opts[locus]; // Increase locus end at 2+4
        }
        *it += offset;
        *it /= (locus_length - 1);
      }
    } 
  }
  
  seg_sites.attr("positions") = wrap(positions); 
  return seg_sites;
}


// [[Rcpp::export]]
List parseSeqgenOutput(const List file_names, 
                       const NumericVector sample_size,
                       const int loci_number,
                       const NumericVector trio_opts = NumericVector(0)) {

  std::string line;
  size_t total_sample_size, locus_length;
  
  List seg_sites(loci_number);
  int locus = -1;
  NumericVector positions(0);
  
  for (int i = 0; i < file_names.size(); ++i) {
    CharacterVector file_name = as<CharacterVector>(file_names(i));
    if (file_name.size() != 1) stop("Expecting one file per locus");
    
    // Open the file
    std::ifstream output(as<std::string>(file_name(0)).c_str(), 
                         std::ifstream::in);
    if (!output.is_open()) {
      stop(std::string("Cannot open file ") + file_name(0));
    }
    
    // We already know the information in the first line
    while ( output.good() ) {
      std::getline(output, line);
      if (line == "") continue;
      if (line.substr(0, 1) == " ") {
        ++locus;
        std::stringstream(line) >> total_sample_size >> locus_length;
        
        seg_sites[locus] = parseSeqgenSegSites(output, locus_length, 
                                               total_sample_size, 
                                               positions, trio_opts);
        
      } else {
        stop("Unexpected line in seq-gen output.");
      }
    }
  
    output.close();
  }
  
  return(seg_sites);
}


