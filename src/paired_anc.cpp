#include <math.h>
#include <fstream>
#include <cmath>



void print_pars(double *par, int n){
  fprintf(stderr, "%f", par[0]);
  for (int i=1;i<n;i++)
    fprintf(stderr, " %f", par[i]);
  fprintf(stderr, "\n");
}

double loglike_paired(double **l, double *pars, int nSites, int nKs){
  double res=0;
  for(int i=0;i<nSites;i++){
    double temp = 0;
    for(int nk=0;nk<nKs;nk++)
      temp += (l[i][nk] * pars[nk]);
    res += log(temp);
  }
  return(res);
}

double prob_gl_anc_af(double *gl1, double f1, double f2){
  double p_g0fa = f1 *f2;
  double p_g1fa = f1 * (1-f2) + f2 * (1-f1);
  double p_g2fa = (1-f1) * (1-f2);
  // double p_g0fa = (1-f1) * (1-f2);
  // double p_g1fa = f1 * (1-f2) + f2 * (1-f1);
  // double p_g2fa = f1 * f2;
  return(gl1[0]*p_g0fa + gl1[1]*p_g1fa + gl1[2]*p_g2fa);
}


void em_anc_paired(double tolStop, int nSites, int nKs, double** pre_calc, double *res2, int &numIter, int nIter){
  int useSq = 1;

  int stepMax = 1;
  int mstep = 4;

  double temppart[nKs];
  double x[nKs];
  double p0[nKs];
  double p1[nKs];
  double q1[nKs];
  double q2[nKs];

  double sr2;
  double sq2;
  double sv2;
  double ttol=0.0000001;
  double rowsum = 0;
  double norma;
  double alpha;


  for(int a=0; a<nKs; a++)
    x[a] = 1.0 / nKs;

  for(int iter=0;iter<nIter;iter++){
    numIter=iter;

    double pars[nKs];
    for(int a=0; a<nKs; a++)
      pars[a] = 0;
    for(int i=0;i<nSites;i++){

      for(int a=0; a<nKs; a++)
        temppart[a] = 0;

      rowsum = 0;

      for(int nk=0;nk<nKs;nk++){
        temppart[nk] = (pre_calc[i][nk] * x[nk]);
        rowsum += temppart[nk];
      }
      for(int nk=0;nk<nKs;nk++)
        pars[nk] += temppart[nk] / rowsum;

    }
    double total = 0;
    for(int nk=0;nk<nKs;nk++)
      total += pars[nk];

    // normalize
    for(int nk=0;nk<nKs;nk++)
      x[nk] = pars[nk] / total;


    // new estimate is now done.

    // if(counter%100==0){
    //   double ll = loglike_paired(pre_calc, x, nSites, nKs);
    //   fprintf(stderr, "%d log: %f ", counter, ll);
    //   print_pars(x, nKs);
    //   fprintf(stderr, "\n");
    // }


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
      for(int j=0;j<nKs;j++)
        p0[j] = x[j];
    }
    if(useSq && iter%3==1){
      for(int j=0;j<nKs;j++)
        p1[j] = x[j];
    }

  }

  for(int j=0;j<nKs;j++)
    res2[j] = x[j];
}

void est_paired_anc(int nSites, int K, int nKs, double *gl1, double **f, double *res2){

  int npop = K;
  double** pre_calc = new double*[nSites];
  for(int i=0;i<nSites;i++){
    pre_calc[i] = new double[nKs];
    int idx = 0;
    for(int a11=0;a11<npop;a11++){
      for(int a12=a11;a12<npop;a12++){
        pre_calc[i][idx] = prob_gl_anc_af(&gl1[i*3], f[i][a11], f[i][a12]);
        // if(i == 0)
        //   fprintf(stderr, "%d %d %f %f %f %f %f %f\n", i, idx, pre_calc[i][idx], f[i][a11], f[i][a12], gl1[i*3+0], gl1[i*3+1], gl1[i*3+2]);
        idx++;
      }
    }
  }




  int maxIter = 5000;
  int currIter = 0;
  double tolStop=0.000001;
  em_anc_paired(tolStop, nSites, nKs, pre_calc, res2, currIter, maxIter);
  double ll = loglike_paired(pre_calc, res2, nSites, nKs);
  fprintf(stderr, "final log (iter: %d): %f ", currIter, ll);
  print_pars(res2, nKs);
  fprintf(stderr, "\n");




  // clean up
  for(int i=0;i<nSites;i++)
    delete[] pre_calc[i];
  delete[] pre_calc;
}