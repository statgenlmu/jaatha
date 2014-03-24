#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

std::vector<double> parseMsPositions(const std::string line);
NumericMatrix parseMsSegSites(std::ifstream &output, const int &seg_sites, const int &individuals);
void addToJsfs(const NumericMatrix seg_sites, const NumericVector sample_size,
               const int &sample_total, NumericMatrix jsfs);

std::string line;

// [[Rcpp::export]]
List parseOutput(const std::string file_name, 
                 const NumericVector sample_size,
                 const int loci_number,
                 const int program = 0,
                 const bool generate_jsfs = true,
                 const bool generate_seg_sites = false,
                 const bool generate_fpc = false) {

  std::ifstream output(file_name.c_str(), std::ifstream::in);
  int individuals = sample_size[0] + sample_size[1];

  if (!output.is_open()) {
    throw exception("Cannot open file");
  }

  NumericMatrix seg_sites(0, 0);
  std::vector<double> positions(0);
  int locus = -1;

  List seg_sites_list(loci_number); 
  NumericMatrix jsfs(sample_size[0]+1, sample_size[1]+1);

  // Output
  while( output.good() ) {
    std::getline(output, line);
    if (line == "//") locus +=1;

    if (line.substr(0, 11) == "segsites: 0") {
      positions.clear(); 
      seg_sites = NumericMatrix(0, 0);
      if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
    } 

    if (line.substr(0, 9) == "segsites:") {
      std::getline(output, line);
      positions = parseMsPositions(line);
      seg_sites = parseMsSegSites(output, positions.size(), individuals);
      if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
      if (generate_jsfs) addToJsfs(seg_sites, sample_size, individuals, jsfs);
    }
  }

  Rprintf("%i\n", positions.size());

  output.close();

  // Return the summary statistics as list without growing it
  List sum_stats;
  if (generate_seg_sites & (!generate_jsfs) & (!generate_fpc)) {
    sum_stats = List::create( _["seg.sites"] = seg_sites_list ) ;
  }
  else if (generate_seg_sites & generate_jsfs & (!generate_fpc)) {
    sum_stats = List::create( _["seg.sites"] = seg_sites_list,
                              _["jsfs"] = jsfs ) ;
  }
  else if ((!generate_seg_sites) & generate_jsfs & (!generate_fpc)) {
    sum_stats = List::create( _["jsfs"] = jsfs ) ;
  }

  return sum_stats;
}

NumericMatrix parseMsSegSites(std::ifstream &output, const int &positions, const int &individuals) {
  NumericMatrix seg_sites(individuals, positions);
  
  for (size_t i = 0; i < individuals; ++i) {
    std::getline(output, line);
    for (size_t j = 0; j < positions; ++j) {
      seg_sites(i,j) = (line[j] == '1');
    }
  }

  return seg_sites;
}

// [[Rcpp::export]]
std::vector<double> parseMsPositions(const std::string line) {
  if (line.substr(0, 11) != "positions: ") {
    throw exception("Failed to read positions from ms' output");
  } 

  std::stringstream stream(line);
  std::vector<double> data;

  // Remove the 'positions: ' at the line's beginning
  stream.ignore(11);

  // Convert the positions into doubles
  std::copy(std::istream_iterator<double>(stream),
            std::istream_iterator<double>(),
            std::back_inserter(data));

  return(data);
}

void addToJsfs(const NumericMatrix seg_sites, const NumericVector sample_size,
               const int &sample_total, 
               NumericMatrix jsfs) {
  int idx1;
  int idx2;
  
  for (int j = 0; j < seg_sites.ncol(); ++j) {
    idx1 = 0;
    idx2 = 0;

    for (int i = 0; i < sample_size[1]; ++i) idx1 += seg_sites(i,j); 
    for (int i = sample_size[1]; i < sample_total; ++i) idx2 += seg_sites(i,j); 

    ++jsfs(idx1, idx2);
  }
}
