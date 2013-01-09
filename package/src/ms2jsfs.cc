/*****************************************************
 *  msFile2jsfs.cc
 *  Calculates the JSFS out of a piped ms output.
 *  
 *  Author:     Lisha Naduvilezhath
 *  Date:       2012-10-05
 *  Licence:    GNU GPLv3 (or later)
 *****************************************************/

#include <R.h>
#include<string>
#include<cstring>
#include<vector>


extern "C"{
  using namespace std;

  //usage: ms2jsfs msoutput samplesize1 samplesize2 anzLoci jsfsArray
  void ms2jsfs(char **msout, int *ps1, int *ps2, int *pnloc, int *jsfs) {
    int s1=*ps1, s2=*ps2, nloc=*pnloc, loc=0, z=5;
    
    //initialize jsfs
    for(int i=0; i<=s1; ++i) {
      for(int j=0; j<=s2; ++j){
	jsfs[j*(s1+1)+i]=0;
      }
    }
    
    //read in locus by locus if there are >0 segsites 
    z=4;
    while (loc < nloc){
      if (strcmp(msout[z],"segsites: 0")){
	int line=z+2;
	int L=strlen(msout[line]); //number of segsites
	for(int p=0; p<L; ++p) {  // go through position by position
	  int x=0, y=0;
	  for(int j=0; j<s1; ++j) {  //pop1
	    x+=(msout[j+line][p]=='1');
	  }
	  for(int j=s1; j<s1+s2; ++j) {  //pop2
	    y+=(msout[j+line][p]=='1');
	  }
	  ++jsfs[y*(s1+1)+x];
	}
	z+=s1+s2+4;
      }
      ++loc;
    }
  }
}
