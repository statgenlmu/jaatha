#include<R.h>	
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<cstring>
#include<vector>
#include<map>
#include<algorithm>

inline unsigned int argmax(int * p, int len) {
  int maxIndex=0;
  for (unsigned int i=1; i<len; ++i) {
    if(p[i]>p[maxIndex]) maxIndex=i;
  }	
  return maxIndex;
};

using namespace std;

/* This program reads in the sequence file that was returned by seqGen and 
   returns the 23 JSFS sumStats:  #positions with tripple polymorphisms, quadrupel 
   polymorphisms, transitions (within and in both pop) and transversions (within and in both pop), 
   and newMig (=bases that are common (>=90%) in one pop and uncommon in the other (<=10%) -> indication for recent migration). */

extern "C" {
  // usage: seqGenFile2jsfs seqgenOutput samplesize1 samplesize2 anzLoci resArraySize resArray
  // 'res' (returned array) will contain the jsfs and possibly additional ss 
  void seqFile2jsfs(char **filename, int *ps1, int *ps2, int *pnloc, int *rSize, int *res){
    //Rprintf("Starting\n");
    ifstream datei(filename[0]);
    //Rprintf("File %d %d\n", filename[0]);
    string currentSeq, seqName_tmp;
    unsigned int numSeq, seqLen, seqName, nLoci=*pnloc, 
                 resSize=*rSize;   
    //Rprintf("resSize %d\n", resSize);
    //Rprintf("nLoci %d\n", nLoci);
    int sampleS[2]={*ps1,*ps2};
    //Rprintf("Sample Sizes %d %d\n", sampleS[0], sampleS[1]);
    char outg;
    numSeq=sampleS[0]+sampleS[1]+1;  // +1for outgroup
    //Rprintf("numSeq %d\n", numSeq);
    vector<vector<char> > sequence(numSeq);
    // boolean vector to note which positions are polymorphic, contain 3 or 4 nucls, 
    //or are transversions          "  fileName:   "<<fileName<<
    //vector<bool> poly2(seqLen), poly3(seqLen),poly4(seqLen), tv(seqLen),ts(seqLen);
    //int res[(sampleS[0]+1)*(sampleS[1]+1)];   //+1 bc 0 needs to be included

    //Rprintf("resSize\n");

    //initialize res
    for(int i=0; i<=resSize; ++i){
      res[i]=0;	     		
    }

    if (datei){
      for (int c=0; c<nLoci; c++){
        //Rprintf("New Locus\n");
        //cout<<"**** locus "<<c<<" :"<<endl;
        datei >>  numSeq >> seqLen;  //first 2 numbers of datei are #seq and #positions
        //cout<<"numSeq:   "<<numSeq<<" seqLen: "<<seqLen<<"    ss1: "<<sampleS[0]<<"    ss2: "<<sampleS[1]<<endl;
        sequence.clear();
        // read sequences into 'sequence[][]' (includes outgroup) and simultaneously sort them
        for(int i=0; i<numSeq; ++i) {
          //Rprintf("i %d\n", i);
          datei >> seqName_tmp >> currentSeq;	

          //ms from phyclust adds an "s" to the seqName. So remove it if it is
          //there...
          if (seqName_tmp.compare(1, 1, "s")) seqName_tmp.erase(0,1);
          seqName = atoi(seqName_tmp.c_str());
          //Rprintf("SeqName: %d", seqName);

          //cout<<"seqName:   "<<seqName<<" currentSeq: "<<currentSeq<<endl;
          sequence[seqName-1].resize(seqLen);   // sequences are unsorted in the file!	
          const char* cstr=currentSeq.c_str();
          for(int p=0; p<seqLen; ++p) {
            sequence[seqName-1][p]=cstr[p];
          }
        }

        //look at each position seperately
        for(int pos=0; pos<seqLen; ++pos) {   
          int derivedCount[2]= {0};     // per pop	     
          outg=sequence[numSeq-1][pos];
          //cout<<"**** position "<<pos<<" :"<<outg<<endl;

          // which nucleotides are present at position pos, including outgroup
          for(int i=0; i<sampleS[0]; ++i) {
            if (sequence[i][pos]!=outg){
              ++derivedCount[0];
            }
          }
          for(int i=sampleS[0]; i<(numSeq-1); ++i) {   //without outgroup
            if (sequence[i][pos]!= outg){
              ++derivedCount[1];
            }
          }
          //cout<<"derivedCount:  "<<derivedCount[0]<<"  "<<derivedCount[1]<<endl;

          //jsfs calculation (positions that are same to outgroup will be counted in res[0,0])
          ++res[derivedCount[1]*(sampleS[0]+1)+derivedCount[0]];
        }  // positions  
      }  //nLoci-for
      datei.close(); 	
    } //if file exists
    else{
      cout<< "ERROR: Cannot open file '"<< filename[0]<<"'"<<endl;		  
    }
  }
  /*cout<<"  0 1 2 3 4 5"<<endl;
    for(int i=0; i<=sampleS[0]; ++i){
    cout<<i<<":  ";
    for(int j=0; j<=sampleS[1]; ++j){
    cout<< res[j*(sampleS[0]+1)+i]<<" ";
    }
    cout<<endl;
    }	*/
}    //extern 
