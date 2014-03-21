#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

void parseLine(std::string line, std::vector<int> &jsfs, const int &s1, const int &s2);

// [[Rcpp::export]]
List parseMsOutput(const std::string file_name, 
                   const NumericVector sample_size,
                   const bool jsfs = true,
                   const bool fpc = false) {

  std::ifstream ms_output(file_name.c_str(), std::ifstream::in);

  
  if (!ms_output.is_open()) {
    throw exception("Cannot open file");
  }

  //while ( ms_output.good() ) {
    //getline(myfile, line);
    //parseLine(line, jsfs, sample_size[0], sample_size[1]);
  //}   

  ms_output.close();

  List sum_stats(2);
  sum_stats[0] = 1;
  sum_stats[1] = "BLUB";
  return sum_stats;
}

//readNextSegSites(std::ifstream
