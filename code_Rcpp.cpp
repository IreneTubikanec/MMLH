
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

//------------------------------------------------------------------------------
// Authors: Anna Melnykova, Irene Tubikanec
// Date: 2023-08-18
//
// Description: auxiliary functions required for the implementation
//              of Algorithm MMLH, proposed in the paper:
//
//              Granger Causal Inference in Multivariate Hawkes Processes
//              by Minimum Message Length, by K. Hlavackova-Schindler, A. Melnykova and I. Tubikanec
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Auxiliary search function
//------------------------------------------------------------------------------
// [[Rcpp::export]]
int which_cpp(double val_time, NumericVector vec_times)
{
  int index_k=0;
  int vec_size=vec_times.size();
  while (vec_times[index_k] < val_time && index_k<(vec_size))
  {
     index_k=index_k+1;
  }
  
  return index_k;
};

//------------------------------------------------------------------------------
// Implementation of expressions Aij(t) in Section 4.3 for node i
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix A_i_cpp(NumericVector beta, List y, int i_in)
{
  int p=beta.size();
  int i=i_in-1;
  NumericVector yi=y[i];
  int ni=yi.size();
  NumericMatrix ret(ni, p);
  
  int k_index;
  double sum1;
  NumericVector yj;
  
  for(int j = 0; j < p; j++)
  {
    yj=y[j];
    for(int l=0; l < ni; l++)
    {
      k_index=which_cpp(yi[l],yj);
      sum1=0.0;
      for(int k=0; k<k_index; k++)
      {
        sum1=sum1+exp(-beta[i]*(yi[l]-yj[k]));  
      }
      ret(l,j)=sum1;
    }
  }
  return ret;
};

//------------------------------------------------------------------------------
// Implementation of 2nd derivatives required for the Hessian matrix in Section 4.3
//------------------------------------------------------------------------------
//#H_mu_cpp
// [[Rcpp::export]]
double H_mu_cpp(double mu, NumericVector alpha, NumericVector beta, NumericVector yi, NumericMatrix Ai)
{
  int p=beta.size();
  double sum1;
  double sum2;

  sum1=0.0;
  for(int l = 0; l < yi.size(); l++)
  {
    sum2=0.0;
    for(int j = 0; j < p; j++)
    {
      sum2=sum2+alpha[j]*Ai(l,j);
    }
    sum1=sum1+1/((mu+sum2)*(mu+sum2));
  }
  return sum1;
};

//#H_mu_alphak_cpp
// [[Rcpp::export]]
double H_mu_alphak_cpp(double mu, NumericVector alpha, NumericVector beta, NumericVector yi, int k_in, NumericMatrix Ai)
{
  int p=beta.size();
  double sum1;
  double sum2;
  int k=k_in-1;
  
  sum1=0.0;
  for(int l = 0; l < yi.size(); l++)
  {
    sum2=0.0;
    for(int j = 0; j < p; j++)
    {
      sum2=sum2+alpha[j]*Ai(l,j);
    }
    sum1=sum1+((Ai(l,k))/((mu+sum2)*(mu+sum2)));
  }
  return sum1;
};

//#H_alphak_alpham_cpp
// [[Rcpp::export]]
double H_alphak_alpham_cpp(double mu, NumericVector alpha, NumericVector beta, NumericVector yi, int k_in, int m_in, NumericMatrix Ai)
{
  int p=beta.size();
  double sum1;
  double sum2;
  int k=k_in-1;
  int m=m_in-1;
  
  sum1=0.0;
  for(int l = 0; l < yi.size(); l++)
  {
    sum2=0.0;
    for(int j = 0; j < p; j++)
    {
      sum2=sum2+alpha[j]*Ai(l,j);
    }
    sum1=sum1+((Ai(l,k)*Ai(l,m))/((mu+sum2)*(mu+sum2)));
  }
  return sum1;
};

//------------------------------------------------------------------------------
// Implementation of the Hessian matrix for node i in Section 4.3
//------------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix Hessian_exact_i_cpp(double mu, NumericVector alpha, NumericVector beta, NumericVector yi, int i, NumericMatrix Ai)
{
  int p=beta.size();
  NumericMatrix ret(p+1,p+1);
  
  ret(0,0)=H_mu_cpp(mu,alpha,beta,yi,Ai);
  
  for(int k = 0; k < p; k++)
  {
    ret(0,k+1)=H_mu_alphak_cpp(mu, alpha, beta, yi, k+1, Ai);
    ret(k+1,0)=ret(0,k+1);
    for(int m = 0; m <= k; m++)
    {
      ret(k+1,m+1)=H_alphak_alpham_cpp(mu, alpha, beta, yi, k+1, m+1, Ai);
      if(m!=k)
      {
        ret(m+1,k+1)=ret(k+1,m+1);
      }
    }
  }
  return ret;
};


//------------------------------------------------------------------------------
// Implementation of the negative log-likelihood for node i in Section 4.2
//------------------------------------------------------------------------------
// [[Rcpp::export]]
double Likelihood_i_cpp(double mu, NumericVector alpha, NumericVector beta, List y, int i_in, NumericMatrix Ai)
{
  double sum1;
  double sum2;
  double sum3;
  double sum4;
  int i=i_in-1;
  int p=beta.size();
  NumericVector yi=y[i];
  int ni=yi.size();
  NumericVector yj;
  int nj;
  double ret;
  
  //Part I
  sum1=0.0;
  for(int j = 0; j < p; j++)
  {
    yj=y[j];
    nj=yj.size();
    sum2=0.0;
    for(int k = 0; k < nj; k++)
    {
      sum2=sum2+(1.0-exp(-beta[i]*(yi[ni-1]-yj[k])));
    }
    sum1=sum1+(alpha[j]/beta[i])*sum2;
  }
  
  //Part II
  sum3=0.0;
  for(int l = 0; l < ni; l++)
  {
    sum4=0.0;
    for(int j = 0; j < p; j++)
    {
      sum4=sum4+alpha[j]*Ai(l,j);
    }
    sum3=sum3+log( mu + sum4 );
  }
  
  ret=mu*yi[ni-1]+sum1-sum3;
  
  return ret;
};



