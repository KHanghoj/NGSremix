#include <iostream>
#include <fstream>
#include "Cinterface.h"
#include "paired_anc.h"
#include <string>
#include <iomanip>//used for setting number of decimals in posterior
#include <cstdlib>
#include <zlib.h>
#include <sstream>
#include <cstring>
#include <vector>
#include <algorithm> // std::find
#include <iostream> // stringstream
#include <sys/stat.h>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif
#include <pthread.h>
#include <assert.h>

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif





#include "filereader_and_conversions.h"
#include "extractors.h"
#include "asort.h"
#include "relateAdmix.h"
#include "ibAdmix.h"

using namespace std;

pthread_t *threads = NULL;
pthread_t *threads1 = NULL;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
int NumJobs;
int jobs;
int *printArray;
FILE *fp;
FILE *fp_paired;
int cunt;
int doInbreeding =0;
int do_both_anc = 0;
int MYSEED = 999999;
int numBoot = 0;
int blockSize = 1000;

eachPars *allPars = NULL;

bool usePlink = false;
bool useBeagle = false;
bool useGlf = false;
std::vector<int> sample_keep;

bool COOL_PA = true;

int VERBOSE = 1;

void readDoubleGZ(double **d,int x,int y,const char*fname,int neg){
  fprintf(stdout,"\t-> Opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cant open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==gzgets(fp,buf,lens)){
 	fprintf(stderr,"Error: Only %d sites in frequency file (Check for MAF filter in ngsremix/admixture\n",i);
	exit(0);
    }
    if(neg)
      d[i][0] = 1-atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = 1-atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  if(NULL!=gzgets(fp,buf,lens)){
    fprintf(stderr,"Error: Too many sites in frequency file. Only %d sites in data file\n",x);
    exit(0);

  }
  gzclose(fp);
}

void readDouble(double **d,int x,int y,const char*fname,int neg){
  fprintf(stdout,"\t-> Opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==fgets(buf,lens,fp)){
    	fprintf(stderr,"Error: Only %d individuals in admxiture (-qname) file\n",i);
      exit(0);
    }
    if(neg)
      d[i][0] = 1-atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = 1-atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
   if(NULL!=fgets(buf,lens,fp)){
    fprintf(stderr,"Error: Too many individuals in admixture (-qname) file. Only %d individuals in data file\n",x);
    exit(0);
  }

  fclose(fp); 
}

int readRow(gzFile gz, char *buf, std::string &row){

  while (gzgets(gz, buf, sizeof(buf)) != Z_NULL){
      std::string temp = buf;
      row += temp;
      if(row[row.size()-1] == '\n'){
        row[row.size()-1] = '\0';
        return row.size();
      }
    }
  return row.size();
}

int ncols(std::string & row){
  std::string word ;
  std::istringstream ss(row) ;
  int ret=0 ;
  while (ss >> word) {
    ret++;
  }
  return ret;
}

void init_param(double * start, const int n){
  double sum=0;
  for (int j=0; j<n; j++){
    start[j] = drand48();
    sum += start[j];
  }

  for (int j=0; j<n; j++)
    start[j] /= sum;
}

void readBeagle(const char *fname, myPars *pars){
  const char *delims = "\t \n";
  
  bool hasHeader = true;
  std::vector<std::string> alldata;
  alldata.reserve(100000);
  
  int lens = 4096;
  // int lens = 1000000;
  std::string row;
  row.reserve(lens);
  int nlines=0 ;
  char buf[lens];

  gzFile fp = gzopen(fname, "rb");
  if (fp==Z_NULL){

    fprintf(stderr,"\n\nERROR: '%s' cannot open file: %s\n\n", __FUNCTION__,fname);
    exit(0);
  };

  fprintf(stdout, "\t-> Beagle - Reading from: %s\n",fname);

  while(readRow(fp, buf, row)!=0){
    if(nlines==0 && hasHeader){
      row.clear();
      hasHeader=false;
      continue;
    }
    alldata.push_back(row);
    row.clear();
    nlines++;

    if(VERBOSE){
      if(nlines % 100000 == 0)
        fprintf(stderr, "\t-> Beagle - %d sites processed\r", nlines);
    }
    
  }
  
  int nSites = nlines;
  int nInd = ncols(alldata[0])/3 - 1;

  fprintf(stdout, "\t-> Beagle - %d sites and %d nInd processed\n", nSites, nInd);
  fprintf(stdout, "\t-> Beagle - Transpose from %d X %d*3 to %d X %d*3\n", nSites, nInd, nInd, nSites);

  dMatrix *returnMat = allocDoubleMatrix(nInd,nSites*3);  
  char *major = new char[nSites];
  char *minor = new char[nSites];
  char **ids = new char*[nSites];

  
  
  for(int i=0; i<nSites; i++){
    char * t = strdup(alldata[i].c_str());
    // see https://github.com/KHanghoj/code_snippets/blob/8bf16e703f8eab3fda593b0dcf9aa6506ff16950/code/read_beagle.cpp
    ids[i] = strdup(strtok(t,delims)); // pos
    major[i] = strtok(NULL,delims)[0]; // major
    minor[i] = strtok(NULL,delims)[0]; // minor
    for(int j=0; j<nInd; j++)
      for(int jj=0; jj<3; jj++)
        returnMat->matrix[j][i*3 + jj] = atof(strtok(NULL, delims));
  }
  pars->dataGL = returnMat;
  pars->major = major;
  pars->minor = minor;
  pars->ids = ids;
  pars->nInd = nInd;
  pars->nSites = nSites;
  
  fflush(stdout);  
}


int getK(const char*fname){
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  if(NULL==gzgets(fp, buf,lens)){
    fprintf(stderr,"Increase buffer\n");
    exit(0);
  }
  strtok(buf,delims);
  int K=1;
  while(strtok(NULL,delims))
    K++;
  gzclose(fp);

  return(K);
}


double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


void info(){
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-beagle    [str] name of the gzipped beagle file (for running from genotype likelihoods)\n");
  fprintf(stderr,"\t-plink     [str] prefix of the binary plink file (excluding the .bed; alternative to -beagle for running from called genotypes))\n");
  fprintf(stderr,"\t-fname     [str] path to ancestral population frequencies file (ending in .fopt.gz if estimated with NGSadmix, .P if estimated with ADMIXUTRE)\n"); 
  fprintf(stderr,"\t-qname     [str] path to admixture proportions file (ending in .qopt if estimated with NGSadmix NGSadmix, .Q if estimated with ADMIXTURE)\n"); 
  fprintf(stderr,"\t-o         [str] name of the output file\n\n"); 
  
  fprintf(stderr,"Paired and Parental Ancestry Estimation\n");
  fprintf(stderr,"\t-bothanc 1 [int] Estimates paired and parental ancestries for each individual. Note Relatedness will NOT be estimated. [Default -bothanc 0]\n\n");
  fprintf(stderr,"\t-boot 0 [int] Number of bootstraps [Default 0]\n\n");
  fprintf(stderr,"\t-block 1000 [int] Number of SNPs in block [Default 0]\n\n");
  
  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-P         [int] Number of threads\n");
  fprintf(stderr,"\t-seed      [uint]\n");
  fprintf(stderr,"\t-tol       [float] Lower tolerance for excluding admixture/paired ancestry estimates [Default: 0.001]\n");  
  fprintf(stderr,"\t-select    [Comma separated (1,2,3) and/or range with dash (1-3), 1-based indexes] Individuals to do estimates for\n");
  fprintf(stderr,"\t-notcool 1 [int] Disables paired ancestry. [Default -notcool 0] \n");

  fprintf(stderr,"Legacy:\n"); 
  fprintf(stderr,"\t-F 1       if you want to estimate inbreeding\n"); 
  fprintf(stderr,"\t-autosomeMax 22\t autosome ends with this chromsome\n"); 

  exit(0);

}

void *functionC(void *a) //the a means nothing
{
  int running_job;

  pthread_mutex_lock(&mutex1);
  
  while (jobs > 0) {
    running_job = jobs--;
    pthread_mutex_unlock(&mutex1);

    ////////////////////////////////////////////// not protected
    int c = NumJobs-running_job;
    eachPars p=allPars[c];
    myPars *pars=p.pars;
    int i=p.ind1;
    int j=p.ind2;
    int numIt=0;
    int myBlockSize=0;
    if( numBoot>0 ){
      if( running_job % (numBoot+1) != 0) // dont boopstrap the first
	myBlockSize = blockSize;
    }
    
    if(usePlink){
      if(COOL_PA)
        relateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q_paired[i],pars->Q_paired[j],p.start,pars->F,pars->tol, COOL_PA);
      else
        relateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol, COOL_PA);
    }else if(useBeagle){
      if(COOL_PA)
        ngsrelateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q_paired[i],pars->Q_paired[j],p.start,pars->F,pars->tol, COOL_PA);
      else
        ngsrelateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol, COOL_PA);        

    }
    p.numIter=numIt;
    p.numI[0]=numIt;
    //////////////////////////////////////////////

    pthread_mutex_lock(&mutex1);

    int d = NumJobs-running_job;
    printArray[d]=1;
    if(d%50==0)
      fprintf(stderr,"\r\trunning i1:%d i2:%d",allPars[d].ind1,allPars[d].ind2);  

    while(cunt<NumJobs){
    
      if(printArray[cunt]==0)
	break;
      if(numBoot==0 )
	fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",allPars[cunt].ind1+1,allPars[cunt].ind2+1,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0]);
      else{
	if(cunt % (numBoot+1) == 0 )
	  fprintf(fp,"org\t%d\t%d\t%f\t%f\t%f\t%d\n",allPars[cunt].ind1+1,allPars[cunt].ind2+1,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0]);
	else
	  fprintf(fp,"boot\t%d\t%d\t%f\t%f\t%f\t%d\n",allPars[cunt].ind1+1,allPars[cunt].ind2+1,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0]);
      }

      
        cunt++;
    }

  }
  pthread_mutex_unlock(&mutex1);

  return NULL;
}

void *functionPaired(void *a) //the a means nothing
{
  int running_job;

  pthread_mutex_lock(&mutex1);
  
  while (jobs > 0) {
    running_job = jobs--;
    pthread_mutex_unlock(&mutex1);

    ////////////////////////////////////////////// not protected
    int c = NumJobs-running_job;
    eachPars p=allPars[c];
    myPars *pars=p.pars;
    int i=p.ind1;
    int numIter=0, numIterPar = 0;
    double llh=0, llhPar=0;
    int do_par = 1;
    int do_not_par= 0;
    int myBlockSize=0;
    int seed = 0;
    if( numBoot>0 ){
      if( running_job % (numBoot+1) != 0){ // dont boopstrap the first
	myBlockSize = blockSize;
	seed=running_job+MYSEED;
      }
      
      
    }

    
    if(useBeagle){
      numIter=    est_paired_anc_gl(pars->nSites, pars->K, pars->dataGL->matrix[i], pars->F, p.pairedPars->pairedParameters, do_not_par, llh , myBlockSize,seed);
      numIterPar= est_paired_anc_gl(pars->nSites, pars->K, pars->dataGL->matrix[i], pars->F, p.pairedPars->parentalParameters, do_par, llhPar , myBlockSize,seed);
      } else if(usePlink){
      numIter=    est_paired_anc_gt(pars->nSites, pars->K, pars->data->matrix[i], pars->F,  p.pairedPars->pairedParameters, do_not_par, llh , myBlockSize,seed);
      numIterPar= est_paired_anc_gt(pars->nSites, pars->K, pars->data->matrix[i], pars->F,  p.pairedPars->parentalParameters, do_par, llhPar, myBlockSize,seed);
        } else {
      
      }
    allPars[c].numIter = numIter;
    allPars[c].numIterPar = numIterPar;
    allPars[c].llh = llh;
    allPars[c].llhPar = llhPar;
    //////////////////////////////////////////////
    pthread_mutex_lock(&mutex1);

    int d = NumJobs-running_job;
    printArray[d]=1;
    if(d%10==0)
      fprintf(stderr,"\r\trunning i1:%d",allPars[d].ind1);  

    while(cunt<NumJobs){
    
      if(printArray[cunt]==0)
        break;

       
      // print stuff to file
      int ind_idx = allPars[cunt].ind1;
      fprintf(fp_paired, "%d", ind_idx+1);
      if(numBoot>0){
	if(cunt % (numBoot+1) == 0 )
	  fprintf(fp_paired, " org");
	else
	  fprintf(fp_paired, " boot");
      }
      int K = pars->K;
      int nKs = (K-1)*K / 2 + K;
      fprintf(fp_paired, " %f", allPars[cunt].pairedPars->pairedParameters[0]);
      for (int ii=1;ii<nKs;ii++){
          fprintf(fp_paired, ",%f", allPars[cunt].pairedPars->pairedParameters[ii]);
      }         
      fprintf(fp_paired, " %d", allPars[cunt].numIter);
      fprintf(fp_paired, " %f", allPars[cunt].llh);
      
      nKs = K*2;
      fprintf(fp_paired, " %f", allPars[cunt].pairedPars->parentalParameters[0]);
      for (int ii=1;ii<nKs;ii++){
	fprintf(fp_paired, ",%f", allPars[cunt].pairedPars->parentalParameters[ii]);
      }         
      fprintf(fp_paired, " %d", allPars[cunt].numIterPar);
      fprintf(fp_paired, " %f", allPars[cunt].llhPar);
      fprintf(fp_paired, "\n");
      
      cunt++;
    }
  }
  pthread_mutex_unlock(&mutex1);
  return NULL;
}




void *functionIBadmix(void *a) //the a means nothing
{
  int running_job;

  pthread_mutex_lock(&mutex1);
  
  while (jobs > 0) {
    running_job = jobs--;
    pthread_mutex_unlock(&mutex1);

    ////////////////////////////////////////////// not protected
    int c = NumJobs-running_job;
    eachPars p=allPars[c];
    myPars *pars=p.pars;
    int i=p.ind1;
    int j=p.ind2;
    int numIt=0;
    double llh=0;
      
    ibAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data,pars->Q,p.start,pars->F,pars->tol,i,llh);

	
    
    p.numIter=numIt;
    p.numI[0]=numIt;
    p.llh=llh;
    p.start[1]=llh;
    //////////////////////////////////////////////

    pthread_mutex_lock(&mutex1);

    int d = NumJobs-running_job;
    printArray[d]=1;
    if(d%50==0)
      fprintf(stderr,"\r\t-> running i1:%d llh:%f",allPars[d].ind1,allPars[d].start[1]);  

    while(cunt<NumJobs){
    
      if(printArray[cunt]==0)
	break;

      fprintf(fp,"%d\t%f\t%f\t%d\n",allPars[cunt].ind1, allPars[cunt].start[0], allPars[cunt].start[1], allPars[cunt].numI[0]);
        cunt++;
    }

  }
  pthread_mutex_unlock(&mutex1);

  return NULL;
}

std::vector<int> parse_select_string(std::string & str){
  std::vector<int> result;
  std::stringstream ss(str);
  while( ss.good() ){
    std::string substr;
    std::getline( ss, substr, ',' );
    if(substr.find("-") != std::string::npos){
      std::vector<int> t;
      std::string substr2;
      std::stringstream ss2(substr);
      
      while(std::getline(ss2, substr2, '-')){
        t.push_back(atoi(substr2.c_str()));
      }
    
      if(t.size()!=2){
        fprintf(stderr, "\nBAD_PARSING!! EXITING\n");
        for(size_t j=0; j<t.size(); j++)
          fprintf(stderr, "%d\n", t[j]);
        exit(0);
      }
      for(int j=t[0]; j<=t[1]; j++)
        result.push_back( j );
      
    } else {
      result.push_back( atoi(substr.c_str()) );
    }
  }
  return result;
}

bool vec_contains(std::vector<int> &vec, int &tok){
  return (std::find(vec.begin(), vec.end(), tok)!=vec.end());
}


void fex(const char* fileName){
  FILE *fp = NULL;
  if((fp=fopen(fileName,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fileName);
    exit(0);
  }
  fclose(fp);
}


int main(int argc, char *argv[]){
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);
  
  // print commandline
  for(int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  
  // if not arguments are given print information
  if(argc==1){
    info();
    return 0;
 }
  
  int useSq=1;
  double tolStop=0.000001;
  int maxIter=5000;
  
  int numIter=4000;
  // testing
  // double tol=0.00001;
  double tol=0.001; // this is the lower boundary for admixture proportions and paired ancestry proportions
  const char *outname = "ngsremix.res";
  int autosomeMax = 23;
  string plink_fam;
  string plink_bim;
  string plink_bed;
  std::string beagle_file, glf_file;
  int nThreads = 1;
  bArray *plinkKeep = NULL; //added in 0.97;
  const char* fname = NULL;
  const char* qname = NULL;
  

  //parse arguments
  int argPos=1;
  while(argPos <argc){
    if (strcmp(argv[argPos],"-o")==0){
      outname  = argv[argPos+1]; 
    }
    else if (strcmp(argv[argPos],"-F")==0)
      doInbreeding = atoi(argv[argPos+1]); 
    else if (strcmp(argv[argPos],"-accel")==0)
      useSq = atoi(argv[argPos+1]); 
    else if (strcmp(argv[argPos],"-a")==0){
      autosomeMax = atoi(argv[argPos+1])+1; 
    }
    else if(strcmp(argv[argPos],"-fname")==0 || strcmp(argv[argPos],"-f")==0) 
      fname=argv[argPos+1]; 
    else if(strcmp(argv[argPos],"-qname")==0 || strcmp(argv[argPos],"-q")==0) 
      qname=argv[argPos+1];
    else if(strcmp(argv[argPos],"-nThreads")==0 || strcmp(argv[argPos],"-P")==0) 
      nThreads=atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-plink")==0){
      usePlink = true;
      std::string p_str =string( argv[argPos+1]);
      if(p_str.length()>4){
	std::string ext = p_str.substr(p_str.length()-4,p_str.length());
	if (!ext.compare(".bed")||!ext.compare(".bim")||!ext.compare(".fam")){
	  std::string front = p_str.substr(0,p_str.length()-4);
	  plink_bim = (front+".bim");
	  plink_fam = (front+".fam");
	  plink_bed = (front+".bed");
	}else{
	  plink_bim = (p_str+".bim");
	  plink_fam = (p_str+".fam");
	  plink_bed = (p_str+".bed");	
	}}else{
	plink_bim = (p_str+".bim");
	plink_fam = (p_str+".fam");
	plink_bed = (p_str+".bed");	
     	}
    } else if (strcmp(argv[argPos], "-beagle")==0){
      beagle_file = argv[argPos+1];
      useBeagle = true;
    } else if (strcmp(argv[argPos], "-glf")==0){
      glf_file = argv[argPos+1];  // not implemented yet
      useGlf = true;
    } else if (strcmp(argv[argPos], "-notcool") == 0){   // -notcool 1 disables paired ancestry
      if (atoi(argv[argPos+1])==1)
        COOL_PA = false;
      else
        COOL_PA = true;
    } else if(strcmp(argv[argPos], "-select") == 0){
      std::string str = string(argv[argPos+1]);
      sample_keep = parse_select_string(str);
    } else if(strcmp(argv[argPos], "-seed") == 0){
      MYSEED = atoi(argv[argPos+1]);    
    } else if(strcmp(argv[argPos], "-boot") == 0){
      numBoot = atoi(argv[argPos+1]);    
    } else if(strcmp(argv[argPos], "-block") == 0){
      blockSize = atoi(argv[argPos+1]);    
    } else if(strcmp(argv[argPos], "-tol") == 0){
      tol = atof(argv[argPos+1]);
    } else if(strcmp(argv[argPos], "-bothanc")==0){
      do_both_anc = atoi(argv[argPos+1]);;
    } else{
      printf("\nArgument unknown will exit: %s \n",argv[argPos]);
      info();
      return 0;
    }

    argPos+=2;
  }

  if((useBeagle && usePlink) || (!useBeagle && !usePlink)){
    fprintf(stderr, "\n\nERROR - cannot provide both (or none of) '-beagle' and '-plink' file\n\n");
    info();
    exit(0);
  }
  if(numBoot > 0 && nThreads == 1 ){
    fprintf(stderr, "\n\nERROR - boostrap is only implemented which threading. use option -P e.g. -P 10 \n\n");
    info();
    exit(0);
  }
  
  //check if files exits
  if (!COOL_PA)
    fex(qname);
  fex(fname);

  int numInds;
  int nInd;
  int nSites;
  myPars *pars =  new myPars();  


  srand(time(NULL));
  if (MYSEED ==  999999){
    MYSEED=rand();
  }
  fprintf(stdout,"\t-> Seed is: %d\n",MYSEED);
  srand48(MYSEED);
  // now set variables with   fprintf(stderr, "%f", drand48());
  // using drand48();
  
  
  if(useBeagle){
    // read beagle and transpose to nInd * nsites*3
    readBeagle(beagle_file.c_str(), pars);
    
    nSites = pars->nSites;
    nInd = pars->nInd;
  
  }else if (usePlink){
    //////////////////////////////////////////////////
    //read plink data or beagle
    printf("\t-> Will assume these are the plink files:\n\t\tbed: %s\n\t\tbim: %s\n\t\tfam: %s\n",plink_bed.c_str(),plink_bim.c_str(),plink_fam.c_str());
    numInds = numberOfLines(plink_fam.c_str())-1;//the number of individuals is just the number of lines
    fprintf(stdout,"\t-> Plink file contains %d samples\n",numInds);
    plinkKeep = doBimFile(pars,plink_bim.c_str()," \t",autosomeMax);
    fprintf(stdout,"\t-> Plink file contains %d autosomale SNPs\n",plinkKeep->numTrue);
    fprintf(stdout,"\t-> reading genotypes\n");
    fflush(stdout);
    usiMatrix *tmp = bed_to_usiMatrix(plink_bed.c_str(),numInds,plinkKeep->x);
    fprintf(stdout,"\t-> done allocating and reading plink \n");
    fflush(stdout);
    if(tmp->y==plinkKeep->numTrue){ 
      pars->data = tmp;
    }else{
      fprintf(stdout,"\t-> extractOK (%d %d) ",tmp->x,plinkKeep->numTrue);
      fflush(stdout);

      pars->data = extractOK(plinkKeep, tmp);
      killMatrix(tmp);
      fprintf(stdout," - done \n");
      fflush(stdout);
    }
    killArray(plinkKeep);
    if(0){
      fprintf(stdout,"\t-> sorting ");
      fflush(stdout);
      mysort(pars,0);  // why sorting?
      fprintf(stdout," - done \n");
      fflush(stdout);
    }

    if(pars->data->y != pars->chr->x || pars->position->x != pars->data->y){
      printf("Dimension of data input doesn't have compatible dimensions, program will exit\n");
      printf("Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d\n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
      return 0;
    }

    nSites = pars->data->y;
    nInd = pars->data->x;
      
    
  } else {
    // fucked. no file is used.
    exit(0);
  }

  if (!sample_keep.empty()){
    
    for (int val: sample_keep){
        if (val<1 || val>nInd){
        fprintf(stderr, "-select index: %d is out of sample range (1-based): 1-%d - EXITING\n", val, nInd);
        exit(0);
      }
    }    

    if(sample_keep.size()<2 && !do_both_anc){
      fprintf(stderr, "-select must contain at least two indices - EXITING\n");
      exit(0);
    }
    // make index zero-based
    for(size_t i=0; i<sample_keep.size();i++)
      sample_keep[i] -= 1;
  } 

  int *keepSamples = new int[nInd];
  for(int i=0; i<nInd; i++){
    if(!sample_keep.empty() && !vec_contains(sample_keep, i))
      keepSamples[i] = 0;
    else
      keepSamples[i] = 1;
  }
  
  int K=getK(fname);
  fprintf(stdout,"\t-> K=%d\tnSites=%d\tnInd=%d\n",K,nSites,nInd);
  pars->maxIter=maxIter;
  pars->tol=tol;
  pars->tolStop=tolStop;
  pars->K=K;
  pars->nSites=nSites;
  pars->useSq=useSq;
  if(!do_both_anc)
    fp=fopen(outname,"w");
  double **F =allocDouble(nSites,K);
  double **Q;
  if (!COOL_PA)
      Q=allocDouble(nInd,K);
  pars->F=F;
  pars->Q=Q;
  if (!COOL_PA)
    readDouble(Q,nInd,K,qname,0);
  readDoubleGZ(F,nSites,K,fname,1);

  int nKs_parental = K*2;
  int nKs_paired = ((K-1)*K/2+K);
  double **Q_parental = allocDouble(nInd, nKs_parental);
  double **Q_paired = allocDouble(nInd, nKs_paired);
  pars->Q_parental = Q_parental;
  pars->Q_paired = Q_paired;
  // est both ancestry exit

  if(do_both_anc | COOL_PA){
    time_t t_paired=time(NULL);
    std::string outname2 = strdup(outname);
    outname2 += ".anccoef";
    fp_paired = fopen(outname2.c_str(), "w");  
    // add header
    fprintf(fp_paired, "ind paired_est paired_iter paired_ll parental_est parental_iter parental_ll\n");
    fprintf(stdout,"\t-> Estimating ancestry coefficients. Dumping to %s\n", outname2.c_str());


    if (nThreads==1){
      for (int i=0; i<nInd;i++){
        if(i%10==0)
          fprintf(stderr, "\r\t-> %d/%d paired ancestries estimated", i, nInd);
        if(keepSamples[i]==0)
          continue;
        int numIter=0, numIterPar = 0;
        double llh=0, llhPar=0;
        int do_par = 1;
        int do_not_par= 0;
        if(useBeagle){
          numIter = est_paired_anc_gl(pars->nSites, pars->K, pars->dataGL->matrix[i], pars->F, pars->Q_paired[i], do_not_par, llh,0,0);
          numIterPar = est_paired_anc_gl(pars->nSites, pars->K, pars->dataGL->matrix[i], pars->F, pars->Q_parental[i], do_par, llhPar,0,0);
        } else if(usePlink){
          numIter = est_paired_anc_gt(pars->nSites, pars->K, pars->data->matrix[i], pars->F, pars->Q_paired[i], do_not_par, llh,0,0);
          numIterPar = est_paired_anc_gt(pars->nSites, pars->K, pars->data->matrix[i], pars->F, pars->Q_parental[i], do_par, llhPar,0,0);
        } else {
          // should not happen
        }
        fprintf(fp_paired, "%d", i+1);
        fprintf(fp_paired, " %f", pars->Q_paired[i][0]);
        for (int ii=1;ii<nKs_paired;ii++){
          fprintf(fp_paired, ",%f", pars->Q_paired[i][ii]);
         }         
        fprintf(fp_paired, " %d", numIter);
        fprintf(fp_paired, " %f", llh);

        fprintf(fp_paired, " %f", pars->Q_parental[i][0]);
        for (int ii=1;ii<nKs_parental;ii++){
          fprintf(fp_paired, ",%f", pars->Q_parental[i][ii]);
         }         
        fprintf(fp_paired, " %d", numIterPar);
        fprintf(fp_paired, " %f", llhPar);
        fprintf(fp_paired, "\n");
      } // per ind
    }
    else {    // threading
      if(sample_keep.empty()){
	NumJobs = nInd * (numBoot+1);
	jobs =  nInd * (numBoot+1);
      } else{
	int tmp_nInd = sample_keep.size();
	NumJobs = tmp_nInd * (numBoot+1);
	jobs =  tmp_nInd * (numBoot+1);
      }

      allPars = new eachPars[NumJobs];
      int allPars_counter = 0 ;
      for(int i=0;i<nInd;i++){
        if(keepSamples[i]==0)
          continue;
	for(int b=0 ; b < (numBoot+1) ; b++ ){
	  allPars[allPars_counter].pairedPars = new pairedParsStruct;
	  allPars[allPars_counter].pairedPars->parentalParameters = new double[nKs_parental];
	  allPars[allPars_counter].pairedPars->pairedParameters =  new double[nKs_paired];
	  allPars[allPars_counter].ind1=i;
	  allPars[allPars_counter].ind2=0;
	  allPars[allPars_counter].numIter=0;
	  allPars[allPars_counter].numIterPar=0;
	  allPars[allPars_counter].pars=pars;
	  allPars[allPars_counter].llh=0;
	  allPars[allPars_counter].llhPar=0;

	  
	  allPars_counter++;
	}
      }
      
      printArray=new int[NumJobs];
      for(int c=0;c<NumJobs;c++){
	printArray[c]=0;
      }    
      
      pthread_t thread1[nThreads];
      
      for (int i = 0; i < nThreads; i++)
	pthread_create(&thread1[i], NULL, &functionPaired, NULL);
      
      // Wait all threads to finish
      for (int i = 0; i < nThreads; i++)
	pthread_join(thread1[i], NULL);
      
    }     
    fclose(fp_paired);
    fprintf(stderr,"\n");
    fprintf(stdout, "\t-> %d ancestry estimates took %ld sec.\n", NumJobs, time(NULL)-t_paired);
    // free stuff
    if (do_both_anc){
      for(int j = 0; j < nSites; j++) 
        delete[] F[j];
      delete[] F;
      if (!COOL_PA){ 
        for(int i = 0; i < nInd; i++)
            delete [] Q[i];
        delete[] Q;
      }
      if(usePlink)
        killMatrix(pars->data);
      else if (useBeagle)
        killMatrix(pars->dataGL);
      delete[] allPars;
      delete[] keepSamples;
      if (nThreads!=1){
        for (int i=0; i<nInd;i++){
          delete [] pars->Q_paired[i];
          delete [] pars->Q_parental[i];
        }
        delete[] pars->Q_paired;
        delete[] pars->Q_parental;
        delete[] printArray;
      }
      fprintf(stdout, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
      fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
      fprintf(stdout, "\t[ALL done] results have been outputted to %s\n",outname2.c_str());
      return(0);
    }
  }

  ///////// print header
  if(!doInbreeding)
    fprintf(fp,"ind1\tind2\tk0\tk1\tk2\tnIter\n");

  if(doInbreeding){//inbreeding
    fprintf(fp,"ind\tF\tllh\t\tnIter\n");

    
    if(nThreads==1){
      double *start=new double[3];
      for(int i=0;i<nInd;i++){
	fprintf(stderr,"\r\trunning i1:%d",i);
	start[0] <- 0.02;
	double llh=0;
        
	ibAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->data,pars->Q,start,pars->F,tol,i,llh);
	fprintf(fp,"%d\t%f\t%f\t%d\n",i,start[0],llh,numIter);
      }
    }
    else{

      NumJobs = nInd;
      jobs =  nInd;
      cunt = 0;
      allPars = new eachPars[NumJobs];

      printArray=new int[NumJobs];
      for(int c=0;c<NumJobs;c++){
	printArray[c]=0;
    
	double *start=new double[1];
	int *numI=new int[1];
	start[0]=0.01;
	int i=c;
	int j=0;
      
	allPars[c].start=start;
	allPars[c].ind1=i;
	allPars[c].ind2=j;
	allPars[c].numIter=0;
	allPars[c].numI=numI;
	allPars[c].pars=pars;
      
    }
    pthread_t thread1[nThreads];
    
    for (int i = 0; i < nThreads; i++)
      pthread_create(&thread1[i], NULL, &functionIBadmix, NULL);
    
  // Wait all threads to finish
    for (int i = 0; i < nThreads; i++)
      pthread_join(thread1[i], NULL);
    
    }
  }// done with inbreeding
  else if(nThreads==1){// relatedness no threads
    double *start=new double[3];
    
    fprintf(stderr,"\t-> running i1:0 i2:0");  
    
    for(int i=0;i<nInd-1;i++){
      for(int j=i+1;j<nInd;j++){

        if(keepSamples[i]==0)
          continue;
        if(keepSamples[j]==0)
          continue;

        // GO HERE
        init_param(start, 3);
        fprintf(stderr,"\r\t-> running i1:%d i2:%d",i,j);
        if(usePlink){
          if(COOL_PA)
            relateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->data->matrix[i],pars->data->matrix[j],pars->Q_paired[i],pars->Q_paired[j],start,pars->F,tol, COOL_PA);
          else
            relateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],start,pars->F,tol, COOL_PA);
        }else if(useBeagle){
          if(COOL_PA)
            ngsrelateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q_paired[i],pars->Q_paired[j],start,pars->F,tol, COOL_PA);
          else
            ngsrelateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q[i],pars->Q[j],start,pars->F,tol, COOL_PA);

        }
        fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",i+1,j+1,start[0],start[1],start[2],numIter);
      }
    }
    delete[] start;
  }else{ // with threads (the cool way)

    if(sample_keep.empty()){
      NumJobs = nInd*(nInd-1)/2;// * (numBoot+1);
      jobs =  nInd*(nInd-1)/2;// * (numBoot+1); 
    } else{
      int tmp_nInd = sample_keep.size();
      NumJobs = tmp_nInd*(tmp_nInd-1)/2;// * (numBoot+1);
      jobs =  tmp_nInd*(tmp_nInd-1)/2;//* (numBoot+1);
    }
    cunt = 0;
    allPars = new eachPars[NumJobs];
    int *indMatrix = new int[nInd*(nInd-1)];
    int cunter=0;
    for(int i=0;i<nInd-1;i++){
      for(int j=i+1;j<nInd;j++){
	//	for(int b=0 ; b < (numBoot+1) ; b++ ){
	  if(keepSamples[i]==0)
	    continue;
	  if(keepSamples[j]==0)
	    continue;
	  indMatrix[cunter*2]=i;
	  indMatrix[cunter*2+1]=j;
	  cunter++;
	  //}
      }
    }

    printArray=new int[NumJobs];
    for(int c=0;c<NumJobs;c++){
      printArray[c]=0;
      double *start=new double[3];
      int *numI=new int[1];
      // GO HERE
      init_param(start, 3);
      int i=indMatrix[c*2];
      int j=indMatrix[c*2+1];
      allPars[c].start=start;
      allPars[c].ind1=i;
      allPars[c].ind2=j;
      allPars[c].numIter=0;
      allPars[c].numI=numI;
      allPars[c].pars=pars;
    }
    pthread_t thread1[nThreads];
    for (int i = 0; i < nThreads; i++)
      pthread_create(&thread1[i], NULL, &functionC, NULL);
    
    // Wait all threads to finish
    for (int i = 0; i < nThreads; i++)
      pthread_join(thread1[i], NULL);
    
    delete[] indMatrix;
  }



  // clean
  fprintf(stderr,"\n");
  for(int j = 0; j < nSites; j++) 
    delete[] F[j];
  delete[] F;
  if (!COOL_PA){ 
    for(int i = 0; i < nInd; i++)
        delete [] Q[i];
    delete[] Q;
  }
  fclose(fp);
  if(usePlink)
    killMatrix(pars->data);
  else if (useBeagle)
    killMatrix(pars->dataGL);
  delete[] allPars;

  if(COOL_PA){
    for (int i=0; i<nInd;i++)
      delete [] pars->Q_paired[i];
    delete[] pars->Q_paired;
  }
  
  delete[] keepSamples;
  
  fprintf(stdout, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stdout, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stdout, "\t[ALL done] results have been outputted to %s\n",outname);

  

  return(0);
}

