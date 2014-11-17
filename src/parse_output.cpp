#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

NumericVector parseMsPositions(const std::string line);

NumericMatrix parseMsSegSites(std::ifstream &output, 
                              const NumericVector positions, 
                              const int individuals);

NumericMatrix parseSeqgenSegSites(std::ifstream &output,
                                  size_t locus_length,
                                  const size_t individuals,
                                  NumericVector &position,
                                  const NumericVector trio_opts);
                       
std::string line;

// [[Rcpp::export]]
List parseOutput(const std::vector<std::string> file_names, 
                 const NumericVector sample_size,
                 const int loci_number,
                 const int program = 0,
                 const NumericVector trio_opts = NumericVector(0)) {
  
  size_t individuals = sample_size[0] + sample_size[1];

  List seg_sites(loci_number);

  NumericVector positions(0);
  int locus = -1;

  for (std::vector<std::string>::const_iterator file_name = file_names.begin();
       file_name != file_names.end(); ++file_name) {
         
    // Open the file
    if (file_name->size() == 0) continue;
    std::ifstream output(file_name->c_str(), std::ifstream::in);
    if (!output.is_open()) {
      stop("Cannot open file");
    }
  
    // ms
    if (program == 0) { 
      while( output.good() ) {
        std::getline(output, line);
        if (line == "//") ++locus;

        else if (line.substr(0, 9) == "segsites:") {
          // Read seg_sites
          if (line.substr(0, 11) == "segsites: 0") {
            NumericMatrix ss = NumericMatrix(0, 0);
            ss.attr("positions") = NumericVector(0);
            seg_sites[locus] = ss;
          } else {
            std::getline(output, line);
            positions = parseMsPositions(line);
            seg_sites[locus] = parseMsSegSites(output, positions, 
                                               individuals);
          }
        }
      }
    }
  
    // seq-gen
    else if (program == 1) {
      // We already know the information in the first line
      size_t total_sample_size, locus_length;
    
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
    }
  
    output.close();
  }
  
  // Return the summary statistics
  return seg_sites;
}

NumericMatrix parseMsSegSites(std::ifstream &output, 
                              const NumericVector positions, 
                              const int individuals) {
  
  NumericMatrix seg_sites(individuals, positions.size());
  seg_sites.attr("positions") = positions;
  
  for (int i = 0; i < individuals; ++i) {
    std::getline(output, line);
    for (int j = 0; j < positions.size(); ++j) {
      seg_sites(i,j) = (line[j] == '1');
    }
  }
  
  return seg_sites;
}


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
  
  NumericMatrix seg_sites(individuals, positions.size());
  
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
    
    seg_sites.attr("positions") = wrap(positions); 
  }
  
  return seg_sites;
}


// [[Rcpp::export]]
NumericVector parseMsPositions(const std::string line) {
  std::stringstream stream(line);
  std::vector<double> data;
  
  if (line.substr(0, 11) != "positions: ") {
    stop("Failed to read positions from ms' output");
  }
  
  // Remove the 'positions: ' at the line's beginning
  stream.ignore(11);
  
  // Convert the positions into doubles
  std::copy(std::istream_iterator<double>(stream),
  std::istream_iterator<double>(),
  std::back_inserter(data));
  
  return(wrap(data));
}
