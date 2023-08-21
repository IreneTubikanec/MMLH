
#-------------------------------------------------------------------------------
# Authors: Anna Melnykova, Irene Tubikanec
# Date: 2023-08-18
#
# Description: auxiliary functions required for the implementation
#              of Algorithm MMLH, proposed in the paper:
#
#              Granger Causal Inference in Multivariate Hawkes Processes
#              by Minimum Message Length, by K. Hlavackova-Schindler, A. Melnykova and I. Tubikanec
#-------------------------------------------------------------------------------
  
#-------------------------------------------
# load file
#------------------------------------------- 

Rcpp::sourceCpp("code_Rcpp.cpp")


#-------------------------------------------------------------------------
# computes the likelihood for node i under a given structure gamma
#-------------------------------------------------------------------------

MLE.wrapper_i <- function(parameters, data, beta, gamma,i){
  p <- length(gamma)
  mu <- parameters[1]
  alpha <- as.numeric(gamma)
  if (sum(gamma)>0){
    alpha[alpha==1] <- parameters[2:length(parameters)]
  } else {
    alpha <- numeric(p)*0
  }
  Ai <- A_i_cpp(beta,data,i)
  Li <- Likelihood_i_cpp(mu, alpha, beta, data,i,Ai)
  return(Li)
}

#-------------------------------------------------------------------------
# computes the likelihood plus exponential prior for node i under a given structure gamma
#-------------------------------------------------------------------------

MLE_and_prior.wrapper_i <- function(parameters, data, beta, gamma,i,c){
  p <- length(gamma)
  mu <- parameters[1]
  alpha <- as.numeric(gamma)
  if (sum(gamma)>0){
    alpha[alpha==1] <- parameters[2:length(parameters)]
  } else {
    alpha <- numeric(p)*0
  }
  Ai <- A_i_cpp(beta,data,i)
  Li <- Likelihood_i_cpp(mu, alpha, beta, data,i,Ai)
  cval <- c
  prior_i <- cval*mu+cval*sum(alpha)-log(cval)*(sum(gamma)+1)
  return(Li+prior_i)
}

#-------------------------------------------------------------------------
# computes the determinant of the Hessian for node i under a given structure gamma
#-------------------------------------------------------------------------

Jtheta.i <- function(parameters, data, beta, gamma,i){
  p = length(gamma)
  mu = parameters[1]
  alpha <- as.numeric(gamma)
  if (sum(gamma)>0){
    alpha[alpha==1] <- parameters[2:length(parameters)]
  } 
  Ai<-A_i_cpp(beta,data,i)
  Hess_i <- Hessian_exact_i_cpp(mu,alpha,beta,data[[i]],i,Ai)
  Jtheta = determinant(Hess_i,logarithm = TRUE)$modulus[[1]]
  return(Jtheta)
}



