#include <cmath>
#include <cstdio>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#include "kmin.h"

//genos[indvidual][sites]
//Q[individual][pop]
//freq[s][pop]
//prod[s[freq,1-freq,freq^2,(1-freq)^2]
typedef struct{
  int nSites;
  int *g;
  double **prod;
}opt;


double llh1(int nSites, double F,int* g,double **prod){
  //  fprintf(stderr,"%F log(F):%f\n",F,log(F));
  
  double llh = 0;
  for(int s=0;s<nSites;s++){
    if(g[s]==0)
      llh += log(prod[s][2]*(1-F)+prod[s][0]*F);
    else if(g[s]==1)
      llh += log(2*prod[s][0]*prod[s][1]*(1-F));
    else if(g[s]==2)
      llh += log(prod[s][3]*(1-F)+prod[s][1]*F);
  }

  //  fprintf(stderr,"res:%f\n",llh);
  return llh;
}

double llh2(double x,void *dats){
  opt *o = (opt*) dats;
  
  return llh1(o->nSites,x,o->g,o->prod);  
}


//simple linear seach
double optim(int nSites,double F,int *g,double **prod){

  
  double min=0.01;
  double max=1-min;
  int nstep = 100;
  double step = (max-min)/nstep;
  //  fprintf(stdout,"min:%f max:%f nstep:%d step:%f\n",min,max,nstep,step);
  double mmax = llh1(nSites,min,g,prod);
  double mmaxid = 0;
  for(double x=min+step;x<max;x+=step){
    double tmp = llh1(nSites,x,g,prod);
    // fprintf(stdout,"%f %f\n",x,tmp);
    if(tmp>mmax){
      mmax = tmp;
      mmaxid = x;
    }
   
  }
  //  fprintf(stdout,"x:%f y:%f\n",mmaxid,mmax);
  return mmaxid;
}

//brent method doesn't work
double optim2(int nSites,int *g,double **prod,double *in){
  opt o ;o.nSites = nSites; o.g = g;o.prod =prod;
  double opti = kmin_brent(llh2,0.001,0.999,(void *) &o,0.001,in);
  return opti;
}

//em method, returns llhvalue of optimized F, F is updated by pointer
double optim3(int nSites,int *g,double **prod,double *in,int maxIter,double tole,int nInformativeSites,double *diff,int *nIter){
  //  fprintf(stderr,"[%s] nsites:%d F:%f maxIter:%d tole:%f nInformativeSites:%d\n",__FUNCTION__,nSites,*in,maxIter,tole,nInformativeSites);
  double F = *in;
  
  double llhOld = log(0);//-Inf
  double newF;
  double llh;
  for(*nIter=0;*nIter<maxIter;(*nIter)++){
    llh =0;
    newF =0;
    for(int s=0;s<nSites;s++){
      double llhSite;
      double Fsite=0;
      if(g[s]==0){
	llhSite = prod[s][2]*(1-F)+prod[s][0]*F;
	Fsite = prod[s][0]*F;
      }else if(g[s]==1)
	llhSite = 2*prod[s][0]*prod[s][1]*(1-F);
      else if(g[s]==2){
	llhSite = prod[s][3]*(1-F)+prod[s][1]*F;
	Fsite = prod[s][1]*F;
      }else{
	//fprintf(stderr,"Problem with range of genotype: %d\n",g[s]);
      }
      llh += log(llhSite);
      //      fprintf(stderr,"Fsite/llhSite:%f\n",Fsite/llhSite);
      newF += Fsite/llhSite;
    }

    F = newF/nInformativeSites;
    *diff = fabs(llh-llhOld);
    //    fprintf(stderr,"res[%d]:F=%f llh:%f diff:%f\n",i,F,llh,llh-llhOld);
    if(*diff<tole){
      //fprintf(stderr,"breaking: diff:%f \n",*diff);
      llhOld = llh;
      break;
    }
    
    llhOld = llh;

  }
  *in = F;
  return llhOld;
}

double ibAdmixEM(int nSites,int *g,double **prod,double *in,int& nIter,double tole,int nInformativeSites,double *diff,int maxIter){
  //  fprintf(stderr,"[%s] nsites:%d F:%f maxIter:%d tole:%f nInformativeSites:%d\n",__FUNCTION__,nSites,*in,maxIter,tole,nInformativeSites);
  double F = in[0];
  
  double llhOld = log(0);//-Inf
  double newF;
  double llh;
  int i=0;
  for(i=0;i<maxIter;i++){
    llh =0;
    newF =0;
    for(int s=0;s<nSites;s++){
      double llhSite;
      double Fsite=0;
      if(g[s]==0){
	llhSite = prod[s][2]*(1-F)+prod[s][0]*F;
	Fsite = prod[s][0]*F;
      }else if(g[s]==1)
	llhSite = 2*prod[s][0]*prod[s][1]*(1-F);
      else if(g[s]==2){
	llhSite = prod[s][3]*(1-F)+prod[s][1]*F;
	Fsite = prod[s][1]*F;
      }else{
	//fprintf(stderr,"Problem with range of genotype: %d\n",g[s]);
      }
      llh += log(llhSite);
      //      fprintf(stderr,"Fsite/llhSite:%f\n",Fsite/llhSite);
      newF += Fsite/llhSite;
    }

    F = newF/nInformativeSites;
    *diff = fabs(llh-llhOld);
    //fprintf(stderr,"res[%d]:F=%f llh:%f diff:%f\n",i,F,llh,llh-llhOld);
    if(*diff<tole){
      //fprintf(stderr,"breaking: diff:%f \n",*diff);
      llhOld = llh;
      break;
    }
    
    llhOld = llh;

  }
  nIter=i;
  //  fprintf(stderr,"nIter:%d %d \n",nIter,maxIter);
  in[0] = F;
  return llhOld;
}

void ibAdmix( double tolStop,int nSites,int K,int maxIter,int useSq,int& numIter,iMatrix *genos, double **Q,double *start,double **f,double tol,int theInd,double &llh){
  //  fprintf(stderr,"\n[%s] nSites:%d K:%d dim.genos(%d,%d)\n",__FUNCTION__,nSites,K,genos->x,genos->y);


  
  double F[genos->x];//start stuff
  for(int i=0;i<genos->x;i++) 
    F[i] = 0.02;
  
  double **prod = new double*[genos->y];  
  for(int s=0;s<genos->y;s++)
    prod[s] = new double[4];

  int i=theInd;
    //precalculate prod matrix
  int nInformativeSites =0;
  for(int s=0;s<genos->y;s++){//loop over sites
    prod[s][0] =0;
    for(int k=0;k<K;k++)
      prod[s][0] += Q[i][k] * f[s][k];
    prod[s][0] = 1-prod[s][0];
    prod[s][1] = 1-prod[s][0];
    prod[s][2] = prod[s][0]*prod[s][0];
    prod[s][3] = prod[s][1]*prod[s][1];
    if(genos->matrix[i][s]!=3)
      nInformativeSites++;
  }
  //    double val=llh1(genos->y,F[i],genos->matrix[i],prod);
  //fprintf(stdout,"llh(ind=%d,F=%f)=%f\n",i,F[i],val);
  double diff;
  
  //  double top=optim3(genos->y,genos->matrix[i],prod,&start[0],numIter,tol,nInformativeSites,&diff,&nIter);


  double top;
  start[0] = 0.0001;
  top=ibAdmixEM(genos->y,genos->matrix[i],prod,start,numIter,tol,nInformativeSites,&diff,1);
  
  if(start[0]>0.0001){
    start[0] = 0.01;
    top=ibAdmixEM(genos->y,genos->matrix[i],prod,start,numIter,tol,nInformativeSites,&diff,maxIter);
  }
  
  llh=top;
  //  fprintf(stdout,"%d\t%f\t%f\t%f\t%d\t%d\n",i,start[0],top,diff,numIter,nInformativeSites);
  //fflush(stdout);
    


}
void ibAdmixOld( double tolStop,int nSites,int K,int nIter,int useSq,int& numIter,iMatrix *genos, double **Q,double *start,double **f,double tol,int theInd){
  fprintf(stderr,"\n[%s] nSites:%d K:%d dim.genos(%d,%d)\n",__FUNCTION__,nSites,K,genos->x,genos->y);


  
  double F[genos->x];//start stuff
  for(int i=0;i<genos->x;i++) 
    F[i] = 0.02;
  
  double **prod = new double*[genos->y];  
  for(int s=0;s<genos->y;s++)
    prod[s] = new double[4];

  for(int i=0;i<genos->x;i++){//loop over individuals  
    //precalculate prod matrix
    int nInformativeSites =0;
    for(int s=0;s<genos->y;s++){//loop over sites
      prod[s][0] =0;
      for(int k=0;k<K;k++)
	prod[s][0] += Q[i][k] * f[s][k];
      prod[s][0] = 1-prod[s][0];
      prod[s][1] = 1-prod[s][0];
      prod[s][2] = prod[s][0]*prod[s][0];
      prod[s][3] = prod[s][1]*prod[s][1];
      if(genos->matrix[i][s]!=3)
	nInformativeSites++;
    }
    //    double val=llh1(genos->y,F[i],genos->matrix[i],prod);
    //fprintf(stdout,"llh(ind=%d,F=%f)=%f\n",i,F[i],val);
    double diff;
    int nIter;
    double top=optim3(genos->y,genos->matrix[i],prod,&F[i],numIter,tol,nInformativeSites,&diff,&nIter);
    fprintf(stdout,"%d\t%f\t%f\t%f\t%d\t%d\n",i,F[i],top,diff,nIter,nInformativeSites);
    fflush(stdout);
    
  }
  




}
