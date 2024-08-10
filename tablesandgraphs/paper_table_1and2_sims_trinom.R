# MAP Sim Trinomial

# Load required libraries
library(rjags)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(RBesT)

# Get args from job file: [n_iter, config id]

theargs <- commandArgs(trailingOnly = TRUE)
n_iter <- as.numeric(theargs[1])
cid <- as.numeric(theargs[2])

theme_set(theme_classic())  

# Set seed 
set.seed(1234)


expit <- function(x){
  exp(x)/(1 + exp(x))
}


map_paper_sims <- function(n, lnRRsd = 1, lnRR_mu = log(0.7), p1_true_mu = 0.4, 
                           p1_true_rho = 20, nexp = 10,
                           iter = 1, return_samps = T, return_counts = T){
  
  # if(is.null(logitp1s) & is.null(p1s)){
  #   stop("At least one of logitp1s and p1s must be supplied!")
  # }
  
  
  # # generate probabilities for each experiment: # NOW TAKING PLACE INSIDE ITERATOR
  # rrs <- exp(lnRRs)
  # p11 <- rrs*p1s^2
  # p12 <- p1s - rrs*p1s^2
  # p22 <- 1 - p1s
  
  prior_p1mu_alpha = 1; prior_p1mu_beta = 1; 
  prior_p1rho_alpha = 1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0; 
  prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2;
  
  # nexp <- length(lnRRs)
  # 
  # if(is.null(p1s)){
  #   p1s <- expit(logitp1s)
  # }
  # 
  # if(length(lnRRs) != length(p1s)){stop("Vectors must be of equal length.")}
  
  #lnRRs <- rescale(lnRRs, true_mu = lnRR_mu, newsig = lnRRsd)
  
  # format n, sample size per experiment: 
  if(length(n) == 1 & nexp > 1){
    n <- rep(n, nexp) # make a vector with a value for each
  }else if(length(n) != nexp){
     stop("Length of n must be 1 or equal to nexp.")
  }
  
  # JAGS model
  cat(" 
  model {
   
    # Prior
    mu_p ~ dbeta(prior_p1mu_alpha, prior_p1mu_beta)
    rho_p ~ dgamma(prior_p1rho_alpha, prior_p1rho_beta)
      
    p1_alpha <- mu_p * rho_p
    p1_beta <- rho_p * (1 - mu_p)

    lnRR_mu ~ dnorm(prior_lnRRmu_mu, prior_lnRRmu_tau)
    lnRR_sig ~ dnorm(0, 1) T(0, ) # dunif(prior_lnRRtau_lo, prior_lnRRtau_hi)
    lnRR_tau <- 1/(lnRR_sig*lnRR_sig)
   
    p1plus.new ~ dbeta(p1_alpha, p1_beta)
      
    theta.new ~ dnorm(lnRR_mu, lnRR_tau)
    log(RR.new) <- theta.new

    p11.new <- p1plus.new * p1plus.new * RR.new
    
    # New Study Data #########
    p12.new <- p1plus.new - p11.new
    p22.new <- 1 - p1plus.new
    psnew[1] <- p11.new
    psnew[2] <- p12.new
    psnew[3] <- p22.new
    z.new ~ dmulti(psnew, n.new)

    # Historic Data #########
    for(i in 1:nexp){
      theta[i] ~ dnorm(lnRR_mu, lnRR_tau)
      log(RR[i]) <- theta[i]

      invRR[i] <- 1/RR[i]
      upbound[i] <- min(1, invRR[i]) 
      theta1[i] ~ dbeta(p1_alpha, p1_beta)
      p1plus[i] <- theta1[i] * upbound[i]

      p11[i]<- RR[i]*pow(p1plus[i],2)
	    p12[i]<- p1plus[i] - RR[i]*pow(p1plus[i],2)
	    p22[i]<- 1 - p1plus[i]
	    
	    ps[i,1] <- p11[i]
      ps[i,2] <- p12[i]
      ps[i,3] <- p22[i]
      z[i,1:3] ~ dmulti(ps[i,1:3],n[i])
    }
  }
  
  ", file = "tbl1and2")
  
  # PARAMETERS AND DATA INIT. ##########################
  
  # output sets aside space for following values:
  # 1. lnRRmu_mean: mean lnRR_mu chain. 
  # 2. lnRRmu_CIlen: CI length for lnRR_mu.
  # 3. lnRRmu_covr: coverage of lnRR_mu central mean 
  # 4. lnRRsd_mean: mean of lnRR_sd chain. 
  # 5. lnRRsd_CIlen: CI length for lnRR_sd.
  # 6. lnRRsd_covr: coverage of lnRR_sd true value.
  # 7. mup_mean: mean of p1 mu chain.
  # 8. mup_CIlen: CI length for p1 mu chain. 
  # 9. mup_covr: coverage rate of p1 central mean
  # 10. rhop_mean: mean of p1 rho.
  # 11. rhop_CIlen: CI length for p1 rho chain.
  # 12. rhop_covr: coverage of p1 rho true value. . 
  # 13. RRindivcovr: coverage rate for experiment RR means
  # 14. p1indivcovr: coverage rate for experiment p1 means
  output <- data.frame(matrix(NA, iter, 21))
  colnames(output) <- c("lnRRmu_mean", "lnRRmu_CIlen", "lnRRmu_covr", 
                        "lnRRsd_mean", "lnRRsd_CIlen", "lnRRsd_covr",
                        "mup_mean", "mup_CIlen", "mup_covr", 
                        "rhop_mean", "rhop_CIlen", "rhop_covr",
                        "RRindivcovr", "RRavgpctbias", "p1indivcovr", "p1avgpctbias",
                        "newRR", "newRRCIlen", 
                        "newp1plus", "newp1plusCIlen", "p11hat")
  
  for(i in 1:iter){
    
    # Generate data: 
    lnRRs <- rnorm(nexp, mean = lnRR_mu, sd = lnRRsd)
    rrs <- exp(lnRRs)
    
    p1_alpha_true <- p1_true_mu * p1_true_rho
    p1_beta_true <- p1_true_rho * (1 - p1_true_mu)
    
    p1s <- rbeta(nexp, shape1 = p1_alpha_true, shape2 = p1_beta_true)
    # p1s <- pmin(p1s, 1/rrs) # this seems like it would lead to bunching at 1/rr. Maybe this doesnt matter. 
    p1s <- p1s * pmin(1, 1/rrs) # TRYING THIS OUT: 
    
    # generate probabilities for each experiment: 
    p11 <- rrs*p1s^2
    p12 <- p1s - rrs*p1s^2
    p22 <- 1 - p1s
    
    p <- cbind(p11, p12, p22)
    
    # generate counts: 
    z <- matrix(nrow = nexp, ncol = 3)
    #colnames(z) <- colnames(p)
    for(j in 1:nrow(p)){
      pi <- p[j,]
      
      z[j,] <- as.numeric(rmultinom(1, n, pi))
    }
    
    # Define data for JAGS
    data <- list(n = n, nexp = nexp, z = z, 
                 n.new = n[1], #n1.new = NA, n11.new = NA,
                 prior_p1mu_alpha = prior_p1mu_alpha, prior_p1mu_beta = prior_p1mu_beta, 
                 prior_p1rho_alpha = prior_p1rho_alpha, prior_p1rho_beta = prior_p1rho_beta,
                 prior_lnRRmu_mu = prior_lnRRmu_mu, prior_lnRRmu_tau = prior_lnRRmu_tau,
                 prior_lnRRtau_lo = prior_lnRRtau_lo, prior_lnRRtau_hi = prior_lnRRtau_hi)
    
    # Compile
    jm <- jags.model("tbl1and2", data = data)
    
    # Burnin
    samps <- update(jm, 5000)
    
    # iterations
    samps <- coda.samples(jm, variable.names = c("lnRR_mu", "lnRR_sig", "mu_p", 
                                                 "rho_p", "p1plus", "RR",
                                                 "p1plus.new", "RR.new"),  
                          n.iter = 20000)
    if(i == 1){
      total_samps = samps[[1]]
    }else{
      total_samps = total_samps + samps[[1]]
    }
    
    stuff <- summary(samps)
    
    # core parameter results:
    output[i,]$lnRRmu_mean <- stuff$statistics["lnRR_mu",1]
    output[i,]$lnRRmu_CIlen <- stuff$quantiles["lnRR_mu",5] - stuff$quantiles["lnRR_mu",1]
    output[i,]$lnRRmu_covr <- stuff$quantiles["lnRR_mu",5] > lnRR_mu & stuff$quantiles["lnRR_mu",1] < lnRR_mu
    output[i,]$lnRRsd_mean <- stuff$statistics["lnRR_sig",1]
    output[i,]$lnRRsd_CIlen <- stuff$quantiles["lnRR_sig",5] - stuff$quantiles["lnRR_sig",1]
    output[i,]$lnRRsd_covr <- stuff$quantiles["lnRR_sig",5] > lnRRsd & stuff$quantiles["lnRR_sig",1] < lnRRsd
    output[i,]$mup_mean <- stuff$statistics["mu_p",1]
    output[i,]$mup_CIlen <- stuff$quantiles["mu_p",5] - stuff$quantiles["mu_p",1]
    output[i,]$mup_covr <- stuff$quantiles["mu_p",5] > p1_true_mu & stuff$quantiles["mu_p",1] < p1_true_mu
    output[i,]$rhop_mean <- stuff$statistics["rho_p",1]
    output[i,]$rhop_CIlen <- stuff$quantiles["rho_p",5] - stuff$quantiles["rho_p",1]
    output[i,]$rhop_covr <- stuff$quantiles["rho_p",5] > p1_true_rho & stuff$quantiles["rho_p",1] < p1_true_rho
    
    # average across studies: 
    studyRRrows <- grep("RR\\[", rownames(stuff$statistics), value = TRUE)
    studyp1rows <- grep("p1plus\\[", rownames(stuff$statistics), value = TRUE)
    output[i,]$RRindivcovr <- (mean(stuff$quantiles[studyRRrows,5] > exp(lnRRs) & stuff$quantiles[studyRRrows,1] < exp(lnRRs)))
    output[i,]$p1indivcovr <- (mean(stuff$quantiles[studyp1rows, 5] > p1s & stuff$quantiles[studyp1rows, 1] < p1s))
    output[i,]$RRavgpctbias <- mean((stuff$statistics[studyRRrows, 1] - exp(lnRRs))/exp(lnRRs))
    output[i,]$p1avgpctbias <- mean((stuff$statistics[studyp1rows, 1] - p1s)/p1s)
    
    # data for the produced prior: 
    output[i,]$newRR <- stuff$statistics["RR.new", 1]
    #output[i,]$newRR2point5 <- stuff$quantiles["RR.new", 1]
    #output[i,]$newRR97point5 <- stuff$quantiles["RR.new", 5]
    output[i,]$newRRCIlen <- stuff$quantiles["RR.new", 5] - stuff$quantiles["RR.new", 1]
    output[i,]$newp1plus <- stuff$statistics["p1plus.new", 1]
    output[i,]$newp1plusCIlen <- stuff$quantiles["p1plus.new", 5] - stuff$quantiles["p1plus.new", 1]
    output[i,]$p11hat <- sum(z[,1])/sum(n)
  }
  
  avg_samps = total_samps / iter
  
  if(return_samps & ! return_counts){
    output <- list(output, avg_samps)
  }else if(!return_samps & return_counts){
    output <- list(output, list("z" = z))
  }else if(return_samps & return_counts){
    output <- list(output, avg_samps, list("z" = z))
  }
  return(output)  
}



# Test parameters ########################################################

# iter = 1; nexp = 10; n = 100; p1_alpha = 15; p1_beta = 30; lnRR_mu = 0 
# lnRR_sig = 0.1; prior_p1mu_alpha = 1; prior_p1mu_beta = 1
# prior_p1rho_alpha = .1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0
# prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2


# Get parameters depending on cid ############################

RRsds <- c(0.075, 0.15)
#RRsds <- c(0.3)
Ns <- c(30, 100)
nexps <- c(5, 25)
lnRRmus <- c(log(0.7), log(1.2))

params <- expand.grid(RRsds, Ns, nexps, lnRRmus)
sdev <- params[cid,1]
N <- params[cid,2]
nexp <- params[cid,3]

# overall log RR: 
lnRRmu <- params[cid,4] # log(0.7)

ess.logRR <- expand.grid(RRsds, Ns)
colnames(ess.logRR) <- c('sd', 'n')
ess.logRR$ess <- NA

t0 <- Sys.time()


# RUN SIMS ###############################################

out <- map_paper_sims(n = N, lnRRsd = sdev, lnRR_mu = lnRRmu, return_samps = F, iter = n_iter, 
                     nexp = nexp, return_counts = F, p1_true_mu = 0.4, p1_true_rho = 20)

print(Sys.time() - t0)
Sys.time() - t0

# OUTPUT #####################################################

out$nexp <- nexp
out$N <- N
out$sd <- sdev
out$lnRRmu <- lnRRmu

# ESS: record effective sample size for RR distributions:
# logRRnew <- log(nextdat$RR)
# 
# fit <- automixfit(logRRnew, type = "norm")
# ess.logRR[ess.logRR[,'n']==N & ess.logRR[,'sd']==sd,'ess'] <- ess(fit, sigma = 1)


write.csv(out, paste0("tbl1and2sims_", as.character(cid), ".csv"))

