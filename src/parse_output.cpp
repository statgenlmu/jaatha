#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

NumericVector parseMsPositions(const std::string line);

NumericMatrix parseMsSegSites(std::ifstream &output, 
                              const NumericVector positions, 
                              const int &individuals);

void addToJsfs(const NumericMatrix &seg_sites,
               const NumericVector &sample_size,
               NumericMatrix &jsfs);

void addToFpc(const NumericMatrix &seg_sites, 
              const NumericVector &positions, 
              const NumericVector &breaks_near,
              const NumericVector &breaks_far,
              NumericMatrix &fpc);

std::string line;

// [[Rcpp::export]]
List parseOutput(const std::string file_name, 
                 const NumericVector sample_size,
                 const int loci_number,
                 const int program = 0,
                 const bool generate_jsfs = true,
                 const bool generate_seg_sites = false,
                 const bool generate_fpc = false,
                 const NumericVector fpc_breaks_near = NumericVector(0),
                 const NumericVector fpc_breaks_far = NumericVector(0)) {
  
  if ((!generate_seg_sites) && (!generate_jsfs) && (!generate_fpc)) return List::create();

  //Rprintf("File: %s \n", file_name.c_str());
  std::ifstream output(file_name.c_str(), std::ifstream::in);
  size_t individuals = sample_size[0] + sample_size[1];

  if (!output.is_open()) {
    throw exception("Cannot open file");
  }

  NumericMatrix seg_sites(0, 0);
  NumericVector positions(0);
  int locus = -1;

  if (generate_fpc) {
    if (fpc_breaks_far.size() == 0 || fpc_breaks_near.size() == 0) 
      throw exception("No breaks for fpc sum stats given");
  }

  List seg_sites_list(loci_number); 
  NumericMatrix jsfs(sample_size[0]+1, sample_size[1]+1);
  NumericMatrix fpc(fpc_breaks_near.size()+2, fpc_breaks_far.size()+2);

  // Output
  while( output.good() ) {
    std::getline(output, line);
    if (line == "//") locus += 1;

    if (line.substr(0, 11) == "segsites: 0") {
      seg_sites = NumericMatrix(0, 0);
      if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
      if (generate_fpc) addToFpc(seg_sites, positions, 
                                 fpc_breaks_near, fpc_breaks_far, fpc);
    } 

    else if (line.substr(0, 9) == "segsites:") {
      //Rprintf("Locus %i\n", locus);
      std::getline(output, line);
      positions = parseMsPositions(line);
      seg_sites = parseMsSegSites(output, positions, individuals);
      if (seg_sites.ncol() == 0) throw exception("Failed to parse seg.sites");
      if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
      if (generate_jsfs) addToJsfs(seg_sites, sample_size, jsfs);
      if (generate_fpc) addToFpc(seg_sites, positions, 
                                 fpc_breaks_near, fpc_breaks_far, fpc);
    }
  }

  output.close();

  // Return the summary statistics as list without growing it
  if (generate_seg_sites & (!generate_jsfs) & (!generate_fpc)) {
    return List::create( _["seg.sites"] = seg_sites_list );
  }
  else if ((!generate_seg_sites) & generate_jsfs & (!generate_fpc)) {
    return List::create( _["jsfs"] = jsfs );
  }
  else if ((!generate_seg_sites) & (!generate_jsfs) & (generate_fpc)) {
    return List::create( _["fpc"] = fpc );
  }

  else if (generate_seg_sites & generate_jsfs & (!generate_fpc)) {
    return List::create( _["seg.sites"] = seg_sites_list,
                         _["jsfs"] = jsfs );
  }
  else if (generate_seg_sites & (!generate_jsfs) & (generate_fpc)) {
    return List::create( _["seg.sites"] = seg_sites_list, 
                         _["fpc"] = fpc );
  }
  else if ((!generate_seg_sites) & generate_jsfs & (generate_fpc)) {
    return List::create( _["jsfs"] = jsfs, 
                         _["fpc"] = fpc );
  }

  else {
    return List::create( _["seg.sites"] = seg_sites_list,
                         _["jsfs"] = jsfs, 
                         _["fpc"] = fpc );
  }

}

NumericMatrix parseMsSegSites(std::ifstream &output, 
                              const NumericVector positions, 
                              const int &individuals) {

  NumericMatrix seg_sites(individuals, positions.size());
  List dimnames = List(2);
  dimnames[1] = positions;
  seg_sites.attr("dimnames") = dimnames;
  
  for (int i = 0; i < individuals; ++i) {
    std::getline(output, line);
    for (int j = 0; j < positions.size(); ++j) {
      seg_sites(i,j) = (line[j] == '1');
    }
  }

  return seg_sites;
}

// [[Rcpp::export]]
NumericVector parseMsPositions(const std::string line) {
  if (line.substr(0, 11) != "positions: ") {
    Rprintf("%s\n", line.c_str());
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

  return(wrap(data));
}

