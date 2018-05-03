/*

    
  anders albrecht@binf.ku.dk made this.

*/
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "abcSum.h"
using namespace std;


void abcSum::printArg(FILE *argFile){
  // Adds to the "args" file in the output?
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSum\t%d\n",doSum);
  fprintf(argFile,"\t-doCounts\t%d\tMust choose -doCount 1\n",doCount);
  fprintf(argFile,"\tUse reference allele  (requires -ref)\n");
  fprintf(argFile,"\tMax depth set to\n", maxDepth);
}


void abcSum::getOptions(argStruct *arguments){
  doSum=angsd::getArg("-doSum",doSum,arguments); // Returns -999, 0 or 1. Sets the action for doSum if the argument was supplied
  // in the commandline.

  if(doSum==0)
    return;

  doCount=angsd::getArg("-doCounts",doCount,arguments); // Same as the last one
  char *ref = NULL;
  char *anc = NULL;
  ref=angsd::getArg("-ref",ref,arguments); //
  anc=angsd::getArg("-anc",anc,arguments);
  maxDepth=angsd::getArg("-maxDepth",maxDepth,arguments); // 
  free(ref);
  free(anc);
  //Error messages
  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"Error: bam or soap input needed for -doSum \n");
    exit(0);
  }
  if( doCount==0){
    fprintf(stderr,"Error: -doSums needs allele counts (use -doCounts 1)\n");
    exit(0);
  }

}


// Constructor
abcSum::abcSum(const char *outfiles,argStruct *arguments,int inputtype){
  //Defaults
  doSum=0;
  doCount=0;
  maxDepth =500;
  outfileG=NULL;

  // If no argument then printArg
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSum")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doSum==0){
    shouldRun[index] =0;
    return;
  }
  
  //make output files
  const char* postfix;
  postfix=".lDepth.gz";
  outfileL = aio::openFileGZ(outfiles,postfix);

  postfix=".sDepth.gz";
  outfileS = aio::openFileGZ(outfiles,postfix);

  postfix= ".qsDep.gz";
  outfileQSDep = aio::openFileGZ(outfiles,postfix);

  postfix= ".qsPosi.gz";
  outfileQSPosi = aio::openFileGZ(outfiles,postfix);

  postfix= ".qsIsop.gz";
  outfileQSIsop = aio::openFileGZ(outfiles,postfix);

  postfix= ".qsMapq.gz";
  outfileQSMap = aio::openFileGZ(outfiles,postfix);

  postfix= ".qsBase.gz";
  outfileQSBase = aio::openFileGZ(outfiles,postfix);

  const char* postfix2;
  postfix2=".gDepth";
  outfileG =  aio::openFile(outfiles,postfix2);


  nInd = arguments->nInd;
  // Initialize total depth values to zero for 500 sites
  totalDepth = new suint [maxDepth];
  for(int j=0;j<maxDepth;j++)
    totalDepth[j] =0;

  // Initialize total counts for all bases found at each depth to zero (depth default=500)
  totalCount = new uint64_t[maxDepth];
  for(int j=0;j<maxDepth;j++)
    totalCount[j] =0;

 // Initialize counts for read quality scores (0 to 50) for all possible bases
  qualityByBase = new uint64_t*[5];
  for(int j=0;j<5;j++)
    qualityByBase[j] =new uint64_t[100];

// Initializing counts for quality scores found (0 to 50) at each depth (default=500)
 qualityByDep = new uint64_t*[maxDepth];
  for(int j=0;j<maxDepth;j++)
    qualityByDep[j] =new uint64_t[100];

// Initializing counts for each base at each depth as a pointer list with a size 500 by default by 4 (one for each base)
  countByDepth = new uint64_t*[maxDepth]; 
  for(int j=0;j<maxDepth;j++)
    countByDepth[j] = new uint64_t[4];

  countByDepthRatio = new double*[maxDepth]; 
  for(int j=0;j<maxDepth;j++)
    countByDepthRatio[j] = new double[4];

// Initializing counts for quality scores/mapq found at possible read positions up to 500 in length 

  qualityByPosi = new uint64_t*[500]; 
  for(int j=0;j<500;j++)
    qualityByPosi[j] = new uint64_t[100];

  qualityByIsop = new uint64_t*[500]; 
  for(int j=0;j<500;j++)
    qualityByIsop[j] = new uint64_t[100];

 qualityByMapq = new uint64_t*[maxDepth]; 
  for(int j=0;j<maxDepth;j++)
    qualityByMapq[j] = new uint64_t[100];


// Initializing the values to 0
  for(int j=0;j<maxDepth;j++){

    for(int b=0;b<4;b++){
    countByDepth[j][b]=0;
    countByDepthRatio[j][b]=0.0;
    }

    for(int b=0;b<100;b++){
      qualityByDep[j][b]=0;
        }
  }
  
  for(int j=0;j<500;j++){
    for(int b=0;b<100;b++){
    qualityByPosi[j][b]=0;
    }
  }

  for(int j=0;j<500;j++){
    for(int b=0;b<100;b++){
    qualityByIsop[j][b]=0;
    }
  }

  for(int j=0;j<maxDepth;j++){
    for(int b=0;b<100;b++){
    qualityByMapq[j][b]=0;
    }
  }


  for(int j=0;j<5;j++){
    for(int b=0;b<100;b++){
    qualityByBase[j][b]=0;
    }
  }
}

abcSum::~abcSum(){
    // Print the depth to file and then close it

  if(doSum==0)
    return; 
  for(int j=0;j<maxDepth;j++){
	fprintf(outfileG,"%d\t",totalDepth[j]);
  }

// Avoid divide by zero by setting 0 values to 1 (Then it is 0 depth divided by 1 count  which is still 0)      
for(int j=0;j< maxDepth;j++){
  if (totalCount[j]==0.0){
    totalCount[j]=1.0;
  }

  for( int b = 0; b < 4; b++ ){
    countByDepthRatio[j][b] = (long double)countByDepth[j][b]/(long double)totalCount[j];
      }
  
    gzprintf(outfileS,"%llu\t%llu\t%llu\t%llu\t%llu\t%f\t%f\t%f\t%f\n",j , 
    countByDepth[j][0],countByDepth[j][1],
    countByDepth[j][2],countByDepth[j][3],
    countByDepthRatio[j][0],countByDepthRatio[j][1],
    countByDepthRatio[j][2],countByDepthRatio[j][3]);
}


for(int j=0;j<500;j++){
  // first column is the range of read positions
    for(int k=0;k<100;k++){
 // all other columns correspond to the count found for base quality k from 0 to 60
	  gzprintf(outfileQSPosi,"%llu\t%llu\n",j, qualityByPosi[j][k]);
    }

 }

 for(int j=0;j<500;j++){
  // first column is the range of read positions
    for(int k=0;k<100;k++){
 // all other columns correspond to the count found for base quality k from 0 to 60
	  gzprintf(outfileQSIsop,"%llu\t%llu\n", j, qualityByIsop[j][k]);
    }
 }


 for(int j=0;j< maxDepth;j++){
   // first column is the range of depths
    for(int k=0;k<100;k++){
    // all other columns correspond to the count found for map quality k from 0 to 100
	  gzprintf(outfileQSMap,"%llu\t%llu\n",j, qualityByMapq[j][k]);
    }
    for(int k=0;k<100;k++){
    // all other columns correspond to the count found for base quality k from 0 to 60
	  gzprintf(outfileQSDep,"%llu\t%llu\n",j, qualityByDep[j][k]);
    }

 }

  for(int j=0;j<5;j++){
   // first column is the range of bases
    for(int k=0;k<100;k++){
  // all other columns correspond to the count found for base quality k from 0 to 60
    gzprintf(outfileQSBase,"%c\t%u\n",intToRef[j], qualityByBase[j][k]);
    }
  }



  if(outfileL!=NULL) gzclose(outfileL);
  if(outfileS!=NULL) gzclose(outfileS);
  if(outfileQSPosi!=NULL) gzclose(outfileQSPosi);
  if(outfileQSIsop!=NULL) gzclose(outfileQSIsop);
  if(outfileQSMap!=NULL) gzclose(outfileQSMap);
  if(outfileQSDep!=NULL) gzclose(outfileQSDep);
  if(outfileQSBase!=NULL) gzclose(outfileQSBase);
  if(outfileG!=NULL) fclose(outfileG);



   for(int s=0;s<maxDepth;s++)
    delete[] countByDepth[s];
    delete[] countByDepth;

   for(int s=0;s<maxDepth;s++)
    delete[] countByDepthRatio[s];
  delete[] countByDepthRatio;
/*

  for(int s=0;s<maxDepth;s++)
    delete[] qualityByDep[s];
  delete[] qualityByDep;

  for(int s=0;s<500;s++)
    delete[] qualityByPosi[s];
  delete[] qualityByPosi;

  for(int s=0;s<5;s++)
    delete[] qualityByBase[s];
  delete[] qualityByBase;
*/
  

}


void abcSum::clean(funkyPars *pars){
  // Frees the list-of-lists ('dat') of the chunks stored under extras in funkyPars constructed with fetch() 
  if(doSum==0)
    return; 

  SumStruct *haplo =(SumStruct *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++)
    delete[] haplo->dat[s];
  delete[] haplo->dat;

  delete haplo;


}


void abcSum::printDepth(funkyPars *pars){
  // Prints the depths to file for each position
  SumStruct *haplo =(SumStruct *) pars->extras[index];
  
  for(int s=0;s< pars->numSites;s++) {
    if(pars-> keepSites[s]==0)
      continue;
    
    int n = haplo->depth[s];
    if(n>maxDepth)
      n=maxDepth-1;
    totalDepth[n]++; 

    gzprintf(outfileL,"%s\t%d\t%d\n",header->target_name[pars->refId],pars->posi[s]+1, haplo->depth[s]);
    
  }

  bufstr.l=0;
  
}

void abcSum::print(funkyPars *pars){

  if(doSum==0)
    return;
  
  printDepth(pars);
}

void abcSum::getDepth(funkyPars *pars){
  
  SumStruct *haplo =(SumStruct *) pars->extras[index]; // Takes the SumStruct chunk out of main funkyPars
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;      
    
    int siteCounts[5] = { 0 , 0 , 0 , 0 , 0 };
    for(int i=0;i< pars->nInd;i++){
      
      uint64_t dep=0;
      for( int b = 0; b < 4; b++ ){
	      dep+=pars->counts[s][i*4+b];
      }
      if(dep==0){
	        haplo->dat[s][i]=4; //Haplotype at position s for file i is set to 4 if missing
	        continue;
      }

      for (int k=0; k<dep;k++){
          /* Read position, quality score and actual bases are extracted by current depth */
          int readPosi = pars->chk->nd[s][i]->posi[k];
          int readIsop = pars->chk->nd[s][i]->isop[k];
          int baseQs = pars->chk->nd[s][i]->qs[k];
          int basenum = refToInt[pars->chk->nd[s][i]->seq[k]];

         // The 500x40 array is updated with a count for the current read position and corresponding quality score 
          qualityByPosi[readPosi][baseQs] +=1;
          qualityByIsop[readIsop][baseQs] +=1;
          qualityByBase[basenum][baseQs] += 1;
     
        }
        

      if(dep>=maxDepth){
      dep=maxDepth-1;
      }

      totalCount[dep] +=dep; // the current depth for site s in file i is added to itself
      for( int b = 0; b < 4; b++ ){
        // Add up the counts for each base seperately by their total depth
        countByDepth[dep][b]+=pars->counts[s][i*4+b];
      }

      
      
        for (int k=0; k<dep;k++){
          /* Read position, quality score and actual bases are extracted by current depth */
          int baseQs = pars->chk->nd[s][i]->qs[k];
          int mapQs = pars->chk->nd[s][i]->mapQ[k];
         // The 500x40 array is updated with a count for the current read position and corresponding quality score 
          qualityByDep[dep][baseQs] +=1;  
          qualityByMapq[dep][mapQs] +=1;        
        }
      
      
      if(doSum==1)//random base
	haplo->dat[s][i] = angsd::getRandomCount(pars->counts[s],i,dep);
      else if(doSum==2)//most frequent base, random tie
	haplo->dat[s][i] = angsd::getMaxCount(pars->counts[s],i,dep);
      
      siteCounts[haplo->dat[s][i]]++; // A C T or G is for site s and ALL files
      // is incremented by 1 depending on the chosen haplotype, but it doesn't go anywhere...
    }
    
    //  fprintf(stdout,"%d sfsdfsdf\n",majorminor);      0
    // call major

    int whichMax=0;
    int NnonMis=siteCounts[0];//number of A C G T for a site (non missing)
    for(int b=1;b<4;b++){
      if(siteCounts[b]>siteCounts[whichMax])
	    whichMax=b;
      NnonMis += siteCounts[b];
    }
    
  }
  
}

void abcSum::run(funkyPars *pars){

  if(doSum==0)
    return;

  //allocate haplotype struct
  SumStruct *haplo = new SumStruct;
  haplo->depth= new int[pars->numSites];

  for(int s=0;s<pars->numSites;s++){
    // Adds up the counts for all bases of aligned reads for each site
    haplo->depth[s] = pars->counts[s][0] + pars->counts[s][1] + pars->counts[s][2] + pars->counts[s][3]; //counts comes from -doCount class, Num A, C, G, T for site s
    
  }
  //get reference information to estimate mismatch
  //ref[s]=="A"
  

  haplo->dat = new int*[pars->numSites];  // Sets an array of array of number of sites x number of files 
  for(int s=0;s<pars->numSites;s++)
    haplo->dat[s] = new int[pars->nInd];


  //haploCall struct to pars // add your struct to the chunk // What is index?
  pars->extras[index] = haplo;
  //get haplotypes
  getDepth(pars);



}
