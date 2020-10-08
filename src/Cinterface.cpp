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
int cunt;
int doInbreeding =0;

eachPars *allPars = NULL;

bool usePlink = false;
bool useBeagle = false;
bool useGlf = false;

int VERBOSE = 1;

void *relateWrap(void *a){
  eachPars *p = (eachPars *)a;
  myPars *pars=p->pars;
  int i=p->ind1;
  int j=p->ind2;
  int numIt=0;
  /*
  fprintf(stderr,"%f\n",pars->tolStop);
  fprintf(stderr,"%f\n",pars->tol);
  fprintf(stderr,"%d\n",pars->nSites);
  fprintf(stderr,"%d\n",i);
  fprintf(stderr,"%d\n",j);
  fprintf(stderr,"%d\n",pars->K);
  fprintf(stderr,"%d\n",pars->maxIter);
  fprintf(stderr,"%d\n",pars->useSq);
  fprintf(stderr,"%d\n",pars->data->matrix[i][0]);
 fprintf(stderr,"%d\n",pars->data->matrix[j][0]);
 fprintf(stderr,"%f\n",pars->Q[j][0]);
fprintf(stderr,"%f\n",p->start[0]);
  */
  //if(doInbreeding==0)
  relateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p->start,pars->F,pars->tol);
  //  else
  //    ibAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data,pars->Q,p->start,pars->F,pars->tol,i,pars->likes);
  p->numIter=numIt;
  return NULL;
}

void readDoubleGZ(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"\t-> Opening : %s with x=%d y=%d\n",fname,x,y);
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
 	fprintf(stderr,"Error: Only %d sites in frequency file (maybe increase buffer)\n",i);
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
  fprintf(stderr,"\t-> Opening : %s with x=%d y=%d\n",fname,x,y);
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

    fprintf(stdout,"\n\nERROR: '%s' cannot open file: %s\n\n", __FUNCTION__,fname);
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
      if(nlines % 10000 == 0)
        fprintf(stdout, "\t-> Beagle - %d sites processed\r", nlines);
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
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"can't open:%s\n",fname);
    exit(0);
  }
  int lens=100000 ;
  char buf[lens];
  if(NULL==fgets(buf,lens,fp)){
    fprintf(stderr,"Increase buffer\n");
    exit(0);
  }
  strtok(buf,delims);
  int K=1;
  while(strtok(NULL,delims))
    K++;
  fclose(fp);

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
  fprintf(stderr,"\t-plink name of the binary plink file (excluding the .bed)\n");
  fprintf(stderr,"\t-beagle name of the gzipped beagle file\n");  
  fprintf(stderr,"\t-fname Ancestral population frequencies\n"); 
  fprintf(stderr,"\t-qname Admixture proportions\n"); 
  fprintf(stderr,"\t-o name of the output file\n"); 

  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-P Number of threads\n");
  fprintf(stderr,"\t-F 1\t if you want to estimate inbreeding\n"); 
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
    if(usePlink)
      relateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol);
    else if(useBeagle)
      ngsrelateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol);

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

      fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",allPars[cunt].ind1,allPars[cunt].ind2,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0]);
      //fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\t%d\n",allPars[cunt].ind1,allPars[cunt].ind2,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0],cunt);
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
    //relateAdmix(pars->tolStop,pars->nSites,pars->K,pars->maxIter,pars->useSq,numIt,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],p.start,pars->F,pars->tol);

	
    
    p.numIter=numIt;
    p.numI[0]=numIt;
    p.llh=llh;
    p.start[1]=llh;
    //fprintf(stdout,"## %d\t%f\t%f\t %d\t %d\n",i,p.start[0],p.llh,p.numI[0],p.numIter);
    //////////////////////////////////////////////

    pthread_mutex_lock(&mutex1);

    int d = NumJobs-running_job;
    printArray[d]=1;
    if(d%50==0)
      fprintf(stderr,"\r->running i1:%d llh:%f",allPars[d].ind1,allPars[d].start[1]);  

    while(cunt<NumJobs){
    
      if(printArray[cunt]==0)
	break;

      fprintf(fp,"%d\t%f\t%f\t%d\n",allPars[cunt].ind1, allPars[cunt].start[0], allPars[cunt].start[1], allPars[cunt].numI[0]);
      //fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\t%d\n",allPars[cunt].ind1,allPars[cunt].ind2,allPars[cunt].start[0],allPars[cunt].start[1],allPars[cunt].start[2],allPars[cunt].numI[0],cunt);
        cunt++;
    }

  }
  pthread_mutex_unlock(&mutex1);

  return NULL;
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
  double tol=0.00001;
  const char *outname = "ngsrelateadmix.res";
  int autosomeMax = 23;
  // string geno= "";
  // string pos = "";
  // string chr = "";
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
  
  //check if files exits
  fex(qname);
  fex(fname);

  int numInds;
  int nInd;
  int nSites;
  myPars *pars =  new myPars();  

  if(useBeagle){
    // read beagle and transpose to nInd * nsites*3
    readBeagle(beagle_file.c_str(), pars);
    
    nSites = pars->nSites;
    nInd = pars->nInd;

    // pars->dataGL->print("beagle","test.txt");
    
  }else if (usePlink){
    //////////////////////////////////////////////////
    //read plink data or beagle
    printf("\t-> Will assume these are the plink files:\n\t\tbed: %s\n\t\tbim: %s\n\t\tfam: %s\n",plink_bed.c_str(),plink_bim.c_str(),plink_fam.c_str());
    numInds = numberOfLines(plink_fam.c_str())-1;//the number of individuals is just the number of lines
    plinkKeep = doBimFile(pars,plink_bim.c_str()," \t",autosomeMax);  
    fprintf(stdout,"\t-> Plink file contains %d autosomale SNPs\n",plinkKeep->numTrue);
    fprintf(stdout,"\t-> reading genotypes ");
    fflush(stdout);
    iMatrix *tmp = bed_to_iMatrix(plink_bed.c_str(),numInds,plinkKeep->x);
    fprintf(stdout," - done \n");
    fflush(stdout);
    if(tmp->y==plinkKeep->numTrue){ 
      pars->data = tmp;
    }else{
      fprintf(stdout,"\t-> extractOK (%d %d) ",tmp->x,plinkKeep->numTrue);
      fflush(stdout);
      
      pars->data = extractOK(plinkKeep,tmp);
      killMatrix(tmp);
      fprintf(stdout," - done \n");
      fflush(stdout);
    }
    killArray(plinkKeep);
    fprintf(stdout,"\t-> sorting ");
    fflush(stdout);
    mysort(pars,0);  // why sorting?
    fprintf(stdout," - done \n");
    fflush(stdout);
    // printf("Dimension of genodata:=(%d,%d), positions:=%d, chromosomes:=%d\n",pars->data->x,pars->data->y,pars->position->x,pars->chr->x);
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
  
  int K=getK(qname);
  fprintf(stderr,"\t-> K=%d\tnSites=%d\tnInd=%d\n",K,nSites,nInd);
  pars->maxIter=maxIter;
  pars->tol=tol;
  pars->tolStop=tolStop;
  pars->K=K;
  pars->nSites=nSites;
  pars->useSq=useSq;

  fp=fopen(outname,"w");
  double **F =allocDouble(nSites,K);
  double **Q =allocDouble(nInd,K);
  pars->F=F;
  pars->Q=Q;
  readDouble(Q,nInd,K,qname,0);
  readDoubleGZ(F,nSites,K,fname,1);

  ///////// print header
  /// testing  
  int nKs = ((K-1)*K/2+K);
  double **paired_anc = allocDouble(nInd,nKs);
  for (int i=0; i<nInd;i++){
    est_paired_anc(pars->nSites, K, nKs, pars->dataGL->matrix[i], pars->F, paired_anc[i]);
    fprintf(stderr, "%d", i);
    for (int ii=0;ii<nKs;ii++)
      fprintf(stderr, " %f", paired_anc[i][ii]);
    fprintf(stderr, "\n");
  }

  pars->Q_paired = paired_anc;


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
    fprintf(stderr,"\t->running i1:0 i2:0");  
    
    for(int i=0;i<nInd-1;i++){
      for(int j=i+1;j<nInd;j++){
	start[0]=0.7;
	start[1]=0.2;
        start[2]=0.1;
        fprintf(stderr,"\r->running i1:%d i2:%d",i,j);
        if(usePlink)
          relateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->data->matrix[i],pars->data->matrix[j],pars->Q[i],pars->Q[j],start,pars->F,tol);
        else if(useBeagle)
          ngsrelateAdmix(tolStop,nSites,K,maxIter,useSq,numIter,pars->dataGL->matrix[i],pars->dataGL->matrix[j],pars->Q[i],pars->Q[j],start,pars->F,tol);
        
        fprintf(fp,"%d\t%d\t%f\t%f\t%f\t%d\n",i,j,start[0],start[1],start[2],numIter);
      }
    }
    delete[] start;
  }else{ // with threads (the cool way)

    NumJobs = nInd*(nInd-1)/2;
    jobs =  nInd*(nInd-1)/2;
    cunt = 0;
    allPars = new eachPars[NumJobs];

    int *indMatrix = new int[nInd*(nInd-1)];
    int cunter=0;
    for(int i=0;i<nInd-1;i++){
      for(int j=i+1;j<nInd;j++){
	indMatrix[cunter*2]=i;
	indMatrix[cunter*2+1]=j;
	cunter++;
      }
    }

    printArray=new int[NumJobs];
    for(int c=0;c<NumJobs;c++){
      printArray[c]=0;
    
      double *start=new double[3];
      int *numI=new int[1];
      start[0]=0.7;
      start[1]=0.2;
      start[2]=0.1;
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
  
  for(int i = 0; i < nInd; i++)
    delete [] Q[i];
  delete[] Q;
 
  fclose(fp);
  if(usePlink)
    killMatrix(pars->data);
  else if (useBeagle)
    killMatrix(pars->dataGL);
  delete[] allPars;
 

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr, "\t[ALL done] results have been outputted to %s\n",outname);

  

  return(0);
}

