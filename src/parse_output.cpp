#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

NumericVector parseMsPositions(const std::string line);

NumericMatrix parseMsSegSites(std::ifstream &output, 
                              const NumericVector positions, 
                              const int individuals);
                              
NumericMatrix parseSeqgenSegSites(std::ifstream &output,
                                  const size_t loci_length,
                                  const size_t individuals,
                                  NumericVector &positions);

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
    stop("Cannot open file");
  }

  NumericMatrix seg_sites(0, 0);
  NumericVector positions(0);
  int locus = -1;

  if (generate_fpc) {
    if (fpc_breaks_far.size() == 0 || fpc_breaks_near.size() == 0) 
      stop("No breaks for fpc sum stats given");
  }

  List seg_sites_list(loci_number); 
  NumericMatrix jsfs(sample_size[0]+1, sample_size[1]+1);
  NumericMatrix fpc(fpc_breaks_near.size()+2, fpc_breaks_far.size()+2);

  // ms
  if (program == 0) { 
    while( output.good() ) {
      std::getline(output, line);
      if (line == "//") ++locus;

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
        if (seg_sites.ncol() == 0) Rf_error("Failed to parse seg.sites");
        if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
        if (generate_jsfs) addToJsfs(seg_sites, sample_size, jsfs);
        if (generate_fpc) addToFpc(seg_sites, positions, 
                                   fpc_breaks_near, fpc_breaks_far, fpc);
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
        
        seg_sites = parseSeqgenSegSites(output, locus_length, 
                                        total_sample_size, positions);
        
        if (generate_seg_sites) seg_sites_list[locus] = seg_sites;
        if (generate_jsfs) addToJsfs(seg_sites, sample_size, jsfs);
        if (generate_fpc) addToFpc(seg_sites, positions, 
                                   fpc_breaks_near, fpc_breaks_far, fpc);
      } else {
        stop("Unexpected line in seq-gen output.");
      }
    }
    
    if (locus != loci_number-1) 
      stop("Unexpected number of loci in seq-gen output.");
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
                              const int individuals) {

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


NumericMatrix parseSeqgenSegSites(std::ifstream &output,
                                  const size_t locus_length,
                                  const size_t individuals,
                                  NumericVector &position) {

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
  
    for (std::vector<double>::iterator it = positions.begin(); it != positions.end(); ++it) {
      *it /= (locus_length - 1);
   }
  
    List dimnames = List(2);
    position = wrap(positions);
    dimnames[1] = position;
    seg_sites.attr("dimnames") = dimnames; 
  }
  
  return seg_sites;
}


// [[Rcpp::export]]
NumericVector parseMsPositions(const std::string line) {
  if (line.substr(0, 11) != "positions: ") {
    stop("Failed to read positions from ms' output");
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

