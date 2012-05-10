#include <R.h>	
#include<fstream>
#include<string>
#include<cstring>
#include<vector>


extern "C"{
  using namespace std;
  
  //usage: msFile2jsfs msoutput samplesize1 samplesize2 anzLoci jsfsArray
  void msFile2jsfs(char **filename, int *ps1, int *ps2, int *pnloc, int *jsfs) {
    int s1=*ps1, s2=*ps2, nloc=*pnloc, loc=0, z=4;
    
    //initialize jsfs
    for(int i=0; i<=s1; ++i){
      for(int j=0; j<=s2; ++j){
	jsfs[j*(s1+1)+i]=0;
      }
    }
    
    ifstream datei;   
    datei.open(filename[0]);
    string msout[(nloc)*(4+s1+s2)+2];

    //read in entire file
    int i=0;
    while(i< ((nloc)*(4+s1+s2)+2) & !datei.eof()){
      getline(datei,msout[i]);
      ++i;
    }
    
    //read in locus by locus if there are >0 segsites 
    while (loc < nloc){ 
      if (strcmp(msout[z].c_str(),"segsites: 0")){
	int line=z+2;
	int L=strlen(msout[line].c_str()); //number of segsites
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

    datei.close();
  }
}
