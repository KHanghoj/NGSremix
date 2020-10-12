#include <math.h>
#include <fstream>
#include <cmath>

int is_missing(double *ary){
  if(fabs(ary[0] - ary[1])<1e-6 && fabs(ary[0] - ary[2])<1e-6 && fabs(ary[1] - ary[2])<1e-6)
    return 1;
  else
    return 0;
}

bool is_nan(double x) { return x != x; }

double ** alloc_and_populate_anc_paired(int K, double *a1){

  double ** res = new double*[K];
  for(int a11=0;a11<K;a11++){
    res[a11] = new double[K];
  }
  int idx = 0;
  for(int a11=0;a11<K;a11++){
    for(int a12=a11;a12<K;a12++){
      if(a11 == a12){
        res[a11][a12] = a1[idx];
      } else {
        res[a11][a12] = a1[idx] / 2;
        res[a12][a11] = a1[idx] / 2;
      }
      idx++;
    }
  }
  return res;
}

void dealloc_anc_paired(int K, double **a){
  for(int a11=0;a11<K;a11++)
    delete[] a[a11];
  delete[] a;
}


void relateAdmix(double tolStop,int nSites,int K,int nIter,int useSq,int& numIter,int *geno1,int *geno2, double *a1,double *a2,double *start,double **f,double tol, bool cool){


  int print=0;
  int totSites=0;
  int *keepSites = new int[nSites];
  int npop=K;


  double ** a1_paired; 
  double ** a2_paired;
  if(cool){
    a1_paired = alloc_and_populate_anc_paired(K, a1);
    a2_paired = alloc_and_populate_anc_paired(K, a2);
  }

  double anc_pair_denom[4];

  if(cool){
    for(int a=0; a<4;a++)
      anc_pair_denom[a] = 0;

    for(int z1=0; z1<2; z1++){
      for(int z2=0; z2<2; z2++){
        for(int a11=0;a11<npop;a11++){
          for(int a12=0;a12<npop;a12++){
            for(int a21=0;a21<npop;a21++){
              if(z1==1 && a11!=a21) // if two diff ancestral pops. k1 (ordered) is 0
                continue;
              for(int a22=0;a22<npop;a22++){
                if(z2==1 && a12!=a22) // if two diff ancestral pops. k2 (ordered) is 0
                  continue;
          
                if((a11==a21 || z1 == 0) && (a12==a22 || z2 == 0))
                  anc_pair_denom[z1*2+z2] += a1_paired[a11][a12] * a2_paired[a21][a22];
              }
            }
          }
        }
      }
    }
  }

  
  for(int i=0;i<nSites;i++){
    if(geno1[i]==3 || geno2[i]==3)
      keepSites[i] = 0;
    else{
      keepSites[i] = 1;
      totSites++;
    }
  }
  if(print){
    for(int k=0;k<K;k++){
      fprintf(stderr,"a: K=%d\t%f\t%f\n",k,a1[k],a2[k]);
    }
    fprintf(stderr,"start = %f\t%f\t%f\n",start[0],start[1],start[2]);
    for(int i=0;i<nSites;i++){
      for(int k=0;k<K;k++){
	if(f[i][k]>1 | f[i][k] <0)
	  fprintf(stderr,"freq is fucked\n");
      }
    }
    fprintf(stderr,"\ntotSites=%d\n",totSites);
  }



 
int* g11=new int[nSites];
int* g21=new int[nSites];
int* g12=new int[nSites];
int* g22=new int[nSites];




double* tempPart=new double[totSites*3];
for(int i=0;i<totSites*3;i++)
  tempPart[i] =0;
 


////////// prep for speed 1
double* Pm11=new double[nSites*npop];
double* Pm12=new double[nSites*npop];
double* Pm21=new double[nSites*npop];
double* Pm22=new double[nSites*npop];



 for(int i=0;i<nSites;i++){
   
   // probablity of data
   g11[i]=0;
   g21[i]=0;
   g12[i]=0;
   g22[i]=0;

   if(geno1[i]<1)
     g11[i] = 1;
   if(geno1[i]<2)
     g12[i] = 1;
   if(geno2[i]<1)
     g21[i] = 1;
   if(geno2[i]<2)
     g22[i] = 1;
   for(int a11=0;a11<npop;a11++){
    if(g11[i]==1)
      Pm11[nSites*a11+i] = (g11[i]-f[i][a11]);
    else
      Pm11[nSites*a11+i] = (g11[i]+f[i][a11]);

    if(g21[i]==1)
      Pm12[nSites*a11+i] = (g21[i]-f[i][a11]);
    else
      Pm12[nSites*a11+i] = (g21[i]+f[i][a11]);

    if(g12[i]==1)
      Pm21[nSites*a11+i] = (g12[i]-f[i][a11]);
    else
      Pm21[nSites*a11+i] = (g12[i]+f[i][a11]);

    if(g22[i]==1)
      Pm22[nSites*a11+i] = (g22[i]-f[i][a11]);
    else
      Pm22[nSites*a11+i] = (g22[i]+f[i][a11]);
    //    fprintf(stderr,"fuck %f\t%d\t%f\n",Pm11[a11*nSites+i],g11[i],f[i][a11]);
    //   if(Pm11[a11*nSites+i]>100)
    //  Rprintf("fuck %f\t%d\t%f\n",Pm11[a11*nSites+i],g11[i],f[a11][i]);
   }
}
 // Rprintf("sadfsadf %f\tg=%d\t%f\n",Pm11[npop*1+9],g11[9],f[1][9]);



//////////done prep 1

// return(ans);


 if (!cool){
for(int a11=0;a11<npop;a11++){
  if(a1[a11] <  tol)
    continue;
for(int a12=0;a12<npop;a12++){
  if(a1[a12] <  tol)
    continue;
for(int a21=0;a21<npop;a21++){
  if(a2[a21] <  tol)
    continue;
for(int a22=0;a22<npop;a22++){
  if(a2[a22] <  tol)
    continue;
   double Pa0 = a1[a11]*a1[a12]*a2[a21]*a2[a22];
   double Pa1 = a1[a11]*a1[a12]*a2[a21]*a2[a22];
   double Pa2 = a1[a11]*a1[a12]*a2[a21]*a2[a22];
   double Sum1=0;
   double Sum2=0;
   for(int k=0;k<npop;k++){
      Sum1+=a1[k]*a2[k];
      Sum2+=a1[k]*a2[k];
   }
   Pa1 /=Sum1;
   Pa2  = Pa1/Sum2;
   int k1keep =1;
   int k2keep =1;
   if(a11!=a21)
     k1keep=0;
   if(a12!=a22)
     k2keep=0;
   int count=0;
  for(int i=0;i<nSites;i++){
   // probablity of data
 
    if(keepSites[i]==0)    
      continue;
    
    tempPart[0*totSites + count] += Pa0*Pm11[a11*nSites+i]*Pm12[a21*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];
    if(g11[i]==g21[i] &k1keep)
      tempPart[1*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=1 k2=0
    //tempPart[2*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=1 k2=0
    if(g12[i]==g22[i] & k2keep){
      tempPart[1*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm12[a21*nSites+i]*Pm21[a12*nSites+i];//k1=0 k2=1
      if(g11[i]==g21[i] &k1keep)
	tempPart[2*totSites + count] += Pa2*Pm11[a11*nSites+i]*Pm21[a12*nSites+i];//k1=1 k2=1
    }
    count++;

  }
 }}}}
 } else if (cool){
for(int a11=0;a11<npop;a11++){
for(int a12=0;a12<npop;a12++){
  if(a1_paired[a11][a12] <  tol) // || a1[a12] > 1-tol)
    continue;
for(int a21=0;a21<npop;a21++){
for(int a22=0;a22<npop;a22++){
  if(a2_paired[a21][a22] <  tol) // || a1[a22] > 1-tol)
    continue;

  double Pa00 = a1_paired[a11][a12] * a2_paired[a21][a22] / anc_pair_denom[0*2+0];
  double Pa01 = a1_paired[a11][a12] * a2_paired[a21][a22] / anc_pair_denom[0*2+1];
  double Pa10 = a1_paired[a11][a12] * a2_paired[a21][a22] / anc_pair_denom[1*2+0];
  double Pa11 = a1_paired[a11][a12] * a2_paired[a21][a22] / anc_pair_denom[1*2+1];
    
  
  int k1keep =1;
  int k2keep =1;
  if(a11!=a21)
    k1keep=0;
  if(a12!=a22)
    k2keep=0;
  int count=0;
  for(int i=0;i<nSites;i++){
    // probablity of data
    
    if(keepSites[i]==0)    
      continue;
    
    tempPart[0*totSites + count] += Pa00*Pm11[a11*nSites+i]*Pm12[a21*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];
    if(g11[i] == g21[i] && k1keep){
      tempPart[1*totSites + count] += Pa10*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=1 k2=0
    }

    if(g12[i] == g22[i] && k2keep){
      tempPart[1*totSites + count] += Pa01*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=0 k2=1

      if(g11[i]==g21[i] && k1keep)
	tempPart[2*totSites + count] += Pa11*Pm11[a11*nSites+i]*Pm21[a12*nSites+i];//k1=1 k2=1     
      
    }
    // tempPart[0*totSites + count] += Pa0*Pm11[a11*nSites+i]*Pm12[a21*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];    
    // if(g11[i]==g21[i] &k1keep)
    //   tempPart[1*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=1 k2=0
    // //tempPart[2*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm21[a12*nSites+i]*Pm22[a22*nSites+i];//k1=1 k2=0
    // if(g12[i]==g22[i] & k2keep){
    //   tempPart[1*totSites + count] += Pa1*Pm11[a11*nSites+i]*Pm12[a21*nSites+i]*Pm21[a12*nSites+i];//k1=0 k2=1
    //   if(g11[i]==g21[i] &k1keep)
    //     tempPart[2*totSites + count] += Pa2*Pm11[a11*nSites+i]*Pm21[a12*nSites+i];//k1=1 k2=1
    // }
    count++;
    
  } // sites
  
 }}}}   // anc pops
 } // cool
  
 

 
 int count = 0;
 for(int i=0;i<nSites;i++){
   int mult = 1;
   if(geno1[i] == 1)
     mult *=2;
   if(geno2[i] == 1)
     mult *=2;
  
     if(keepSites[i]==0)    
       continue;

     tempPart[0*totSites+count] *=mult;
     tempPart[1*totSites+count] *=mult;
     tempPart[2*totSites+count] *=mult;
  
     if(geno1[i] == 1 && geno2[i] == 1){
       tempPart[1*totSites+count] /=2;
       tempPart[2*totSites+count] /=2; 
     }
     count++;

 }
 //return(ans);
//// the em part


 int stepMax = 1;
 int mstep = 4;
 int stepMin = 1;
 
 double x[3];
 

 double p0[3];
 double p1[3];
 double q1[3];
 double q2[3];
 double sr2;
 double sq2;
 double sv2;
 double ttol=0.0000001;
 double norma;
 double alpha;
 double siteSum;
 double Pr0;
 double Pr1;
 double Pr2;

 ////// Test if k0>0.999
 x[0]=0.999;
 x[1]=0.0005;
 x[2]=0.0005;
 Pr0 = x[0];
 Pr1 = x[1]/2;
 Pr2 = x[2];

 for(int j=0;j<3;j++)
   x[j] = 0;
 
 for(int i=0;i<totSites;i++){
   
   siteSum=tempPart[i]*Pr0 +
     tempPart[1*totSites + i]*Pr1+
     tempPart[2*totSites + i]*Pr2;
   
   x[0] += tempPart[i]*Pr0/siteSum;
   x[1] += tempPart[1*totSites + i]*Pr1/siteSum;
   x[2] += tempPart[2*totSites + i]*Pr2/siteSum;
 }
 
 for(int j=0;j<3;j++)
   x[j] /= totSites;
 
 ///// start accelerated EM if k0<0.999
 if(x[1]<0.0005 && x[2]<0.0005 ){
   tolStop=0;
   numIter=1;
 }
 else{
   for(int j=0;j<3;j++)
     x[j] = start[j];
   
   for(int iter=0;iter<nIter;iter++){
     numIter=iter;

 
     Pr0 = x[0];
     Pr1 = x[1]/2;
     Pr2 = x[2];
     
     for(int j=0;j<3;j++)
       x[j] = 0;
     
     for(int i=0;i<totSites;i++){
       
       siteSum=tempPart[i]*Pr0 +
	 tempPart[1*totSites + i]*Pr1+
	 tempPart[2*totSites + i]*Pr2;
       
       x[0] += tempPart[i]*Pr0/siteSum;
       x[1] += tempPart[1*totSites + i]*Pr1/siteSum;
       x[2] += tempPart[2*totSites + i]*Pr2/siteSum;
     }


     for(int j=0;j<3;j++)
       x[j] /= totSites;
     //////  acceleration
     if(useSq && iter%3==2){
    
       sr2=0;
       sq2=0;
       sv2=0;

       //get stop sizes
       for(int j=0;j<3;j++){
	 q1[j] = p1[j] - p0[j];
	 sr2+= q1[j]*q1[j];
	 q2[j] = x[j] - p1[j];
	 sq2+= q2[j]*q2[j];
	 sv2+= (q2[j]-q1[j])*(q2[j]-q1[j]); 
       }

       //Stop the algorithm if the step size less than tolStop
       if(sqrt(sr2)<tolStop || sqrt(sq2)<tolStop || (p1[0]>0.999 & q2[0]>0)){
	 tolStop=sr2;
	 break;
       }
       
       //calc alpha and map into [1,stepMax] if needed
       alpha = sqrt(sr2/sv2);
       if(alpha>stepMax)
	 alpha=stepMax;
       if(alpha<1)
	 alpha=1;

       //the magical step
       for(int j=0;j<3;j++)
	 x[j] = p0[j] + 2 * alpha * q1[j] + alpha*alpha * (q2[j] - q1[j]);
       
       //in the rare instans that the boundarys are crossed. map into [ttol,1-ttol]
       for(int j=0;j<3;j++){
	 if(x[j]<ttol)
	   x[j]=ttol;
	 if(x[j]>1-ttol)
	   x[j]=1-ttol;

       }
       norma=x[0]+x[1]+x[2];
       for(int j=0;j<3;j++)
	 x[j] /= norma;

       //change step size
       if (alpha == stepMax) 
	 stepMax = mstep * stepMax;
     }
     
     if(useSq && iter%3==0){
       for(int j=0;j<3;j++)
	 p0[j] = x[j];
     }
     if(useSq && iter%3==1){
       for(int j=0;j<3;j++)
	 p1[j] = x[j];
     }
     
   }

 }


for(int j=0;j<3;j++)
  start[j] = x[j];

delete[] tempPart;
delete[] g11;
delete[] g12;
delete[] g21;
delete[] g22;
delete[] Pm11;
delete[] Pm21;
delete[] Pm12;
delete[] Pm22;
delete[] keepSites;

 if(cool){
   dealloc_anc_paired(K, a1_paired);
   dealloc_anc_paired(K, a2_paired); 
 }

 
}


void ngsrelateAdmix(double tolStop,int nSites,int K,int nIter,int useSq,int& numIter, double *gl1, double *gl2, double *a1,double *a2,double *start,double **f,double tol, bool cool){

  
  // make 2d matrix of the Q paired for each indi;
  
  double ** a1_paired; 
  double ** a2_paired;
  if(cool){
    a1_paired = alloc_and_populate_anc_paired(K, a1);
    a2_paired = alloc_and_populate_anc_paired(K, a2);
  }
  // fprintf(stderr, "test\n\n");
  // for(int a11=0;a11<K;a11++){
  //   for(int a12=0;a12<K;a12++){
  //     fprintf(stderr, "%f ", a1_paired[a11][a12]);
  //   }
  //   fprintf(stderr, "\n");
  // }
  // exit(0);

  int print=0;
  int totSites=0;
  int *keepSites = new int[nSites];
  int npop=K;


  double anc_pair_denom[4];

  if(cool){
    for(int a=0; a<4;a++)
      anc_pair_denom[a] = 0;

    for(int z1=0; z1<2; z1++){
      for(int z2=0; z2<2; z2++){
        for(int a11=0;a11<npop;a11++){
          for(int a12=0;a12<npop;a12++){
            for(int a21=0;a21<npop;a21++){
              if(z1==1 && a11!=a21) // if two diff ancestral pops. k1 (ordered) is 0
                continue;
              for(int a22=0;a22<npop;a22++){
                if(z2==1 && a12!=a22) // if two diff ancestral pops. k2 (ordered) is 0
                  continue;
          
                if((a11==a21 || z1 == 0) && (a12==a22 || z2 == 0))
                  anc_pair_denom[z1*2+z2] += a1_paired[a11][a12] * a2_paired[a21][a22];
              }
            }
          }
        }
      }
    }
  }
  // for(int a=0; a<4;a++)
  //   fprintf(stderr, "%f ", anc_pair_denom[a]);
  // fprintf(stderr, "\n");

  
  for(int i=0;i<nSites;i++){
    if(is_missing(&gl1[i*3]) || is_missing(&gl2[i*3]))
      keepSites[i] = 0;
    else{
      keepSites[i] = 1;
      totSites++;
    }
  }
  if(print){
    for(int k=0;k<K;k++){
      fprintf(stderr,"a: K=%d\t%f\t%f\n",k,a1[k],a2[k]);
    }
    fprintf(stderr,"start = %f\t%f\t%f\n",start[0],start[1],start[2]);
    for(int i=0;i<nSites;i++){
      for(int k=0;k<K;k++){
	if(f[i][k]>1 || f[i][k] <0)
	  fprintf(stderr,"freq is fucked\n");
      }
    }
    fprintf(stderr,"\ntotSites=%d\n",totSites);
  }

  double* tempPart=new double[totSites*3];
  for(int i=0;i<totSites*3;i++) // k0, k1, k2
    tempPart[i] =0;

for(int a11=0;a11<npop;a11++){
  if(!cool && a1[a11] <  tol) // || a1[a11] > 1-tol)
    continue;
for(int a12=0;a12<npop;a12++){
  if(cool && a1_paired[a11][a12] <  tol) // || a1[a12] > 1-tol)
    continue;
  if(!cool && a1[a12] < tol)
    continue;
for(int a21=0;a21<npop;a21++){
  if(!cool && a2[a21] <  tol) // || a1[a21] > 1-tol)
    continue;
for(int a22=0;a22<npop;a22++){
  if(cool && a2_paired[a21][a22] <  tol) // || a1[a22] > 1-tol)
    continue;
  if(!cool && a2[a22] < tol)
    continue;


for(int z1=0; z1<2; z1++){
  if(z1==1 && a11!=a21) // if two diff ancestral pops. k1 (ordered) is 0
    continue;
for(int z2=0; z2<2; z2++){
  if(z2==1 && a12!=a22) // if two diff ancestral pops. k2 (ordered) is 0
    continue;

  
  // fprintf(stderr, "%d %d %d %d %d %d %f %f %f %f\n", a11, a12, a21, a22, z1, z2, a1_paired[a11][a12], a2_paired[a21][a22],  a1_paired[a11][a12] * a2_paired[a21][a22], denom);

for(int g11=0; g11<2; g11++){  // integrate the unobs ordered genotypes
for(int g12=0; g12<2; g12++){
for(int g21=0; g21<2; g21++){
  if(z1==1 && g11!=g21)
    continue;
for(int g22=0; g22<2; g22++){
  if(z2==1 && g12!=g22)
    continue;


  // double denom = get_denom_paired_anc(npop, a1_paired, a2_paired, z1, z2);
  // double Pa = a1_paired[a11][a12] * a2_paired[a21][a22] / denom;

  // if cool method
  double Pa = 0;
  if(cool){
    Pa =  a1_paired[a11][a12] * a2_paired[a21][a22] / anc_pair_denom[z1*2+z2];
  }else{
    Pa = a1[a11]*a1[a12]*a2[a21]*a2[a22]; // this has to be fixed
    
    double sum1=0, sum2 = 0;
    for(int k=0;k<npop;k++){
      sum1+=a1[k]*a2[k];
      sum2+=a1[k]*a2[k];      
    }
    
    if(z1==1)
      Pa /= sum1;
    
    if(z2==1)
      Pa /= sum2;    
  }
  
  int count=0;
  for(int i=0;i<nSites;i++){
   // probablity of data
    
    if(keepSites[i]==0)    
      continue;
    // kh:
    double pgl1 = gl1[i*3 + g11+g12];
    double pgl2 = gl2[i*3 + g21+g22];

    double P1, P2;
      
    if(z1==1 && g11==g21)
      P1 = g21+pow(-1, g11)*f[i][a11];
    else
      P1 = (g11+pow(-1, g11)*f[i][a11]) * (g21+pow(-1, g21)*f[i][a21]);

    if(z2==1 && g12==g22)
      P2 = g22+pow(-1, g22)*f[i][a22];
    else
      P2 = (g12+pow(-1, g12)*f[i][a12]) * (g22+pow(-1, g22)*f[i][a22]);

    tempPart[(z1+z2) * totSites + count] += Pa * P1 * P2 * pgl1 * pgl2;
    count++;

  }   // sites
 }}}} // g11,g12,g21,g22
 }}   // z1 and z2
 }}}} // npops

//// the em part

 int stepMax = 1;
 int mstep = 4;
 int stepMin = 1;
 
 double x[3];
 

 double p0[3];
 double p1[3];
 double q1[3];
 double q2[3];
 double sr2;
 double sq2;
 double sv2;
 double ttol=0.0000001;
 double norma;
 double alpha;
 double siteSum;
 double Pr0;
 double Pr1;
 double Pr2;

 ////// Test if k0>0.999
 x[0]=0.999;
 x[1]=0.0005;
 x[2]=0.0005;
 Pr0 = x[0];
 Pr1 = x[1]/2;
 Pr2 = x[2];

 for(int j=0;j<3;j++)
   x[j] = 0;
 
 for(int i=0;i<totSites;i++){
   
   siteSum=tempPart[i]*Pr0 + tempPart[1*totSites + i]*Pr1 + tempPart[2*totSites + i]*Pr2;
   x[0] += tempPart[i]*Pr0/siteSum;
   x[1] += tempPart[1*totSites + i]*Pr1/siteSum;
   x[2] += tempPart[2*totSites + i]*Pr2/siteSum;
 }

 for(int j=0;j<3;j++)
   x[j] /= totSites;
 
 ///// start accelerated EM if k0<0.999
 if(x[1]<0.0005 && x[2]<0.0005 ){
   tolStop=0;
   numIter=1;
 }
 else{
   for(int j=0;j<3;j++)
     x[j] = start[j];
   
   for(int iter=0;iter<nIter;iter++){
     numIter=iter;

 
     Pr0 = x[0];
     Pr1 = x[1]/2;
     Pr2 = x[2];
     
     for(int j=0;j<3;j++)
       x[j] = 0;
     
     for(int i=0;i<totSites;i++){
       
       siteSum=tempPart[i]*Pr0 +
	 tempPart[1*totSites + i]*Pr1+
	 tempPart[2*totSites + i]*Pr2;
       
       x[0] += tempPart[i]*Pr0/siteSum;
       x[1] += tempPart[1*totSites + i]*Pr1/siteSum;
       x[2] += tempPart[2*totSites + i]*Pr2/siteSum;
     }


     for(int j=0;j<3;j++)
       x[j] /= totSites;
     //////  acceleration
     if(useSq && iter%3==2){
    
       sr2=0;
       sq2=0;
       sv2=0;

       //get stop sizes
       for(int j=0;j<3;j++){
	 q1[j] = p1[j] - p0[j];
	 sr2+= q1[j]*q1[j];
	 q2[j] = x[j] - p1[j];
	 sq2+= q2[j]*q2[j];
	 sv2+= (q2[j]-q1[j])*(q2[j]-q1[j]); 
       }

       //Stop the algorithm if the step size less than tolStop
       if(sqrt(sr2)<tolStop || sqrt(sq2)<tolStop || (p1[0]>0.999 & q2[0]>0)){
	 tolStop=sr2;
	 break;
       }
       
       //calc alpha and map into [1,stepMax] if needed
       alpha = sqrt(sr2/sv2);
       if(alpha>stepMax)
	 alpha=stepMax;
       if(alpha<1)
	 alpha=1;

       //the magical step
       for(int j=0;j<3;j++)
	 x[j] = p0[j] + 2 * alpha * q1[j] + alpha*alpha * (q2[j] - q1[j]);
       
       //in the rare instans that the boundarys are crossed. map into [ttol,1-ttol]
       for(int j=0;j<3;j++){
	 if(x[j]<ttol)
	   x[j]=ttol;
	 if(x[j]>1-ttol)
	   x[j]=1-ttol;

       }
       norma=x[0]+x[1]+x[2];
       for(int j=0;j<3;j++)
	 x[j] /= norma;

       //change step size
       if (alpha == stepMax) 
	 stepMax = mstep * stepMax;
     }
     
     if(useSq && iter%3==0){
       for(int j=0;j<3;j++)
	 p0[j] = x[j];
     }
     if(useSq && iter%3==1){
       for(int j=0;j<3;j++)
	 p1[j] = x[j];
     }
     
   }

 }


for(int j=0;j<3;j++)
  start[j] = x[j];

delete[] tempPart;
delete[] keepSites;

 if(cool){
   dealloc_anc_paired(K, a1_paired);
   dealloc_anc_paired(K, a2_paired); 
 }
}




