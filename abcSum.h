#include "abc.h" // superClass where run/clean/print is inherited from


typedef struct{

  int *depth; //pointer
  int **dat; //pointer of pointers
}SumStruct; 

class abcSum:public abc{
private:

  //non optional arguments
  int doSum;
  int doCount;

  //  float hello;
  
  gzFile outfileL, outfileS, outfileQSDep, outfileQSPosi,outfileQSIsop,outfileQSMap, outfileQSBase;
  FILE* outfileG; //normal file

  int nInd;
  suint *totalDepth; //pointers
  uint64_t *baseCount;
  uint64_t **qualityByBase;
  uint64_t **qualityByDep;
  uint64_t **qualityByPosi;
  uint64_t **qualityByIsop;
  uint64_t **qualityByMapq;
  uint64_t **countByDepth;
  double **countByDepthRatio;
  
  uint64_t *totalCount; //suint to contain large counts, but will be converted to double in ratio calculations (52-64 bit vs 64 bit)
  int maxDepth;
  

  //functions 
  void printDepth(funkyPars *pars);
  void getDepth(funkyPars *pars);

  //print buffer
  kstring_t bufstr;

  
public:
  abcSum(const char *outfiles,argStruct *arguments,int inputtype); //constructor // Secure
  ~abcSum();  // destructor // Secure
  void getOptions(argStruct *arguments); // is called from  constructor
  void run(funkyPars  *pars); // chunks //insecure
  void print(funkyPars *pars); //chunks // secure
  void clean(funkyPars *pars); //chunks // insecure
  void printArg(FILE *argFile);// called from getOptions

};


