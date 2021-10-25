void ibAdmix( double tolStop,int nSites,int K,int maxIter,int useSq,int& numIter,usiMatrix *genos, double **Q,double *start,double **f,double tol,int theInd,double &llh) ;




typedef struct  { //pars used acceleration
  double tol;
  double tolStop;
  

//global mess for accell
  double alpha;
  int stepMax ;
  int mstep ;
  int stepMin ;

  int maxIter;
  int useSq;
  int K;
  int nSites;
  double likes;

}EMoptions ;//pars used by relateHMM
 
