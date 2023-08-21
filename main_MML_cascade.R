
#-------------------------------------------------------------------------------
# load required packages and files
#-------------------------------------------------------------------------------

library("parallel")
library(MLmetrics)
library(hawkes)
source(file="aux_Functions.R")

#-------------------------------------------------------------------------------
# choose setting
#-------------------------------------------------------------------------------

N <- 100   # number of trials (over which the respective F1-score is averaged) 

p <- 7   # dimension of the MHP
Time <- 200 # time horizon of the MHP

bval <- 10^5 # value for uniform prior
cval <- 10^-5 # value for exponential prior 

beta <- rep(1,p)  # true beta constants (fixed)
mu_const <- 0.5 # true mu parameters
mu <- rep(mu_const, p) 
alpha_const <- 0.55 # true alpha parameters

# number of cores for parallel computation
ncl <- detectCores()

#-------------------------------------------------------------------------------
# generate N trajectories 
#-------------------------------------------------------------------------------

trajs <- list() # list of the N trajectories
true_gammas <- list() # list of the corresponding N true gamma structures
for (i in 1:N){ 
  
  # true_alpha (here, cascade structure with self-excitation in the first component)
  alpha_vec <- rep(alpha_const,p*p) 
  alpha <- matrix(alpha_vec,byrow=TRUE,nrow=p)  
  for(k in 1:p){
    for(j in 1:p){
      if(k!=(j+1)){
        alpha[k,j]<-0
      }
    }
  }
  alpha[1,1] <- alpha_const
  alpha_vec <- c(t(alpha))
  
  # corresponding true gamma structure
  true_gamma <- as.numeric(alpha_vec!=0) 
  true_gammas[[i]] <- true_gamma
  
  # trajectory, using the true alpha
  trajs[[i]] <- simulateHawkes(mu,alpha,beta,Time)
}

#-------------------------------------------------------------------------------
# structure set Gamma
#-------------------------------------------------------------------------------

# generate all possible binary combinations
gammas_full <- expand.grid(replicate(p, 0:1, simplify = FALSE)) 

# remove option with only zeros (if ki>0 is assumed)
gammas <- gammas_full[2:dim(gammas_full)[1],]

#-------------------------------------------------------------------------------
# based on each of the N observed trajectories, estimate the corresponding structure
#-------------------------------------------------------------------------------

results <- mclapply(trajs, function(traj) {
  
  results_mml_unif <- numeric(dim(gammas)[1])  
  results_mml_exp <- numeric(dim(gammas)[1])  
  pred <- numeric()
  pred_exp <- numeric()
  par0 <- list()

  for (j in 1:p){
      start <- (j-1)*p+1 # to navigate in the vector alpha_vec
      end <- j*p
      k <- 0  
      
      # computation of the MML criterion for each combination
      for (i in 1:dim(gammas)[1]){        
        gam <- as.numeric(gammas[i,])
        if (all(gam == true_gamma[start:end])){
          true_index <- i
        }  
        par0[[i]] <- mu_const + runif(1,-0.4,0.4)  # initialization of the parameter vector 
        if (sum(gam)!=0) {
          par0[[i]] <- c(par0[[i]],alpha_const+runif(sum(gam),-0.15,0.15)) 
        }
        
        # theta_hat under uniform prior (=MLE)
        pre_result <- optim(par0[[i]],MLE.wrapper_i, data = traj, beta = beta, gamma = gam, i = j, control = list(maxit = 50, factr = 1e3), method = "Nelder-Mead")
        
        # theta_hat under exponential prior
        pre_result_prior <- optim(par0[[i]],MLE_and_prior.wrapper_i, data = traj, beta = beta, gamma = gam, i = j, c = cval, control = list(maxit = 50, factr = 1e3), method = "Nelder-Mead")
        
        # number of non-zero (=1) entries in the structure gam
        k <- sum(gam)
        
        # MML criterion under uniform prior
        detJ <- Jtheta.i(pre_result$par, data = traj, beta = beta, gamma = gam, i =j)
        results_mml_unif[i] <- pre_result$val+log(bval)*(k+1)+0.5*detJ-k/2*log(2*pi)+0.5*log(k*pi)+log(choose(p,k))+log(p+1)+digamma(1) 
        
        # MML criterion under exponential prior
        detJ_exp <- Jtheta.i(pre_result_prior$par, data = traj, beta = beta, gamma = gam, i =j)
        results_mml_exp[i] <- pre_result_prior$val+0.5*detJ_exp-k/2*log(2*pi)+0.5*log(k*pi)+log(choose(p,k))+log(p+1)+digamma(1) 
    }
  
    # find optimal structure and return it  
    good_index <- which(results_mml_unif == min(results_mml_unif))
    good_index_exp <- which(results_mml_exp == min(results_mml_exp))

    pred[start:end] <- as.numeric(gammas[good_index,])
    pred_exp[start:end] <- as.numeric(gammas[good_index_exp,])
  }
  
   return(list(pred, pred_exp))
}, mc.cores = ncl)

#-------------------------------------------------------------------------------
# Store the results (estimated structures) and the underlying true gammas
#-------------------------------------------------------------------------------

estimates_MML_uniform <- matrix(0,nrow=N,ncol=p^2)
estimates_MML_exp <- matrix(0,nrow=N,ncol=p^2)

ground_truth_gammas <- matrix(0,nrow=N,ncol=p^2)

for(i in 1:N){
  estimates_MML_uniform[i,] <- results[[i]][[1]]
  estimates_MML_exp[i,] <- results[[i]][[2]]

  ground_truth_gammas[i,] <- true_gammas[[i]]
}

#------------------------------------------------
# Compare each of the N estimated gammas with the corresponding true_gamma (under the F1-score)
#------------------------------------------------

f1_mml <- rep(0,N)
f1_mml_exp <- rep(0,N)

for(i in 1:N){
  res_MML <- results[[i]][[1]]
  res_MML_exp <- results[[i]][[2]]

  # true gamma
  true_gamma <- true_gammas[[i]]
  
  #F1 score: MML under uniform prior
  res_MML <- factor(as.character(res_MML), levels=c("0","1"))
  val_MML_0 <- table(res_MML,true_gamma)[1,1]
  val_MML_1 <- table(res_MML,true_gamma)[2,2]
  if(val_MML_1==0){
    f1_mml[i] <- 0
  } else {
    f1_mml[i] <- F1_Score(true_gamma, res_MML, positive="1")
  }
  
  #F1 score: MML under exponential prior
  res_MML_exp <- factor(as.character(res_MML_exp), levels=c("0","1"))
  val_MMLexp_0 <- table(res_MML_exp,true_gamma)[1,1]
  val_MMLexp_1 <- table(res_MML_exp,true_gamma)[2,2]
  if(val_MMLexp_1==0){
    f1_mml_exp[i] <- 0
  } else {
    f1_mml_exp[i] <- F1_Score(true_gamma, res_MML_exp,positive="1")
  }
}

#-------------------------------------------------------------------------------
# Compute the F1-scores averaged over the N trials
#-------------------------------------------------------------------------------

f1_mml_r <- mean(f1_mml)
f1_mml_exp_r <- mean(f1_mml_exp)

# Show resulting F1-scores
f1_mml_r # MML under uniform prior
f1_mml_exp_r # MML under exponential prior





