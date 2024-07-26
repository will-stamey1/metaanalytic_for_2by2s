# MAP Sim

# Load required libraries
library(rjags)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(stringr)
library(RBesT)

theme_set(theme_classic())

# Set seed 
set.seed(1234)

# In these simulations, we will remove noise by holding the data 
# constant. For increases in standard deviation, the differences 
# in sampled parameters will be scaled up. 

expit <- function(x){
  exp(x)/(1 + exp(x))
}

rescale <- function(x, true_mu, newsig = 1){
  #mu <- mean(x)
  mu <- true_mu
  sig <- sd(x)
  
  sigmult <- newsig/sig
  
  for(i in 1:length(x)){
    x[i] <- (x[i] - mu) * sigmult + mu
  }
  
  return(x)
}


map_paper_sims <- function(lnRRs, logitp1s = NULL, p1s = NULL, n, lnRRsd = 1, 
                           lnRR_mu = log(0.7), p1_true_mu = 0.4, p1_true_rho = 20,
                           iter = 1, return_samps = T, return_counts = T){
  
  if(is.null(logitp1s) & is.null(p1s)){
    stop("At least one of logitp1s and p1s must be supplied!")
  }
  
  prior_p1mu_alpha = 1; prior_p1mu_beta = 1; 
  prior_p1rho_alpha = .1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0; 
  prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2;
  
  nexp <- length(lnRRs)
  
  if(is.null(p1s)){
    p1s <- expit(logitp1s)
  }
  
  if(length(lnRRs) != length(p1s)){stop("Vectors must be of equal length.")}
  
  lnRRs <- rescale(lnRRs, true_mu = lnRR_mu, newsig = lnRRsd)
  
  # format n, sample size per experiment: 
  if(length(n) == 1 & nexp > 1){
    n <- rep(n, nexp) # make a vector with a value for each
  }else if(length(n) != nexp){
    stop("Length of n must be 1 or equal to nexp.")
  }
  
  #  JAGS model
  cat(" 
  model {
   
    # Prior
    mu_p ~ dbeta(prior_p1mu_alpha, prior_p1mu_beta)
    rho_p ~ dgamma(prior_p1rho_alpha, prior_p1rho_beta)
      
    p1_alpha <- mu_p*rho_p
    p1_beta <- rho_p * (1 - mu_p)
      
    lnRR_mu ~ dnorm(prior_lnRRmu_mu, prior_lnRRmu_tau)
    lnRR_sig ~ dunif(prior_lnRRtau_lo, prior_lnRRtau_hi)
    lnRR_tau <- 1/(lnRR_sig*lnRR_sig)
   
    for(i in 1:nexp){  
      p2[i] <- p1[i] * RR[i] 
    }
    p2.new <- p1.new * RR.new
      
    # New Study Data #########
    n1.new ~ dbin(p1.new, n.new)
    n11.new ~ dbin(p2.new, n1.new)
    
    p1.new ~ dbeta(p1_alpha, p1_beta)
      
    theta.new ~ dnorm(lnRR_mu, lnRR_tau)
    log(RR.new) <- theta.new
    
    # Historic Data #########
    for(i in 1:nexp){
      n1[i] ~ dbin(p1[i], n[i])
      n11[i] ~ dbin(p2[i], n1[i])
      
      theta1[i] ~ dbeta(p1_alpha, p1_beta)
      invRR[i] <- 1/RR[i]
      upbound[i] <- min(1, invRR[i])  
      p1[i] <- theta1[i] * upbound[i]

      log(RR[i]) <- theta[i]
      theta[i] ~ dnorm(lnRR_mu, lnRR_tau)
    }
  }
  
  ", file = "one_samp_zero")
  
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
  output <- data.frame(matrix(NA, iter, 23))
  colnames(output) <- c("lnRRmu_mean", "lnRRmu_CIlen", "lnRRmu_covr", 
                        "lnRRsd_mean", "lnRRsd_CIlen", "lnRRsd_covr",
                        "mup_mean", "mup_CIlen", "mup_covr", 
                        "rhop_mean", "rhop_CIlen", "rhop_covr",
                        "RRindivcovr", "RRavgpctbias", "p1indivcovr", "p1avgpctbias",
                        "newRR", "newRR2point5", "newRR97point5", "newRRCIlen", 
                        "newp1", "newp1CIlen", "pi11hat")
  
  for(i in 1:iter){
    
    # Generate data: 
    n1 <- rbinom(nexp, n, prob = p1s)
    n11 <- rbinom(nexp, n1, prob = (p1s * exp(lnRRs))) 
    
    # Define data for JAGS
    data <- list(n = n, nexp = nexp, n1 = n1, n11 = n11, 
                 n.new = 0, n1.new = 0, n11.new = 0,
                 prior_p1mu_alpha = prior_p1mu_alpha, prior_p1mu_beta = prior_p1mu_beta, 
                 prior_p1rho_alpha = prior_p1rho_alpha, prior_p1rho_beta = prior_p1rho_beta,
                 prior_lnRRmu_mu = prior_lnRRmu_mu, prior_lnRRmu_tau = prior_lnRRmu_tau,
                 prior_lnRRtau_lo = prior_lnRRtau_lo, prior_lnRRtau_hi = prior_lnRRtau_hi)
    
    # Compile
    jm <- jags.model("one_samp_zero", data = data)
    
    # Burnin
    samps <- update(jm, 5000)
    
    # iterations
    samps <- coda.samples(jm, variable.names = c("lnRR_mu", "lnRR_sig", "mu_p", 
                                                 "rho_p", "p1", "RR",
                                                 "p1.new", "RR.new"),  
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
    studyp1rows <- grep("p1\\[", rownames(stuff$statistics), value = TRUE)
    output[i,]$RRindivcovr <- (mean(stuff$quantiles[studyRRrows,5] > exp(lnRRs) & stuff$quantiles[studyRRrows,1] < exp(lnRRs)))
    output[i,]$p1indivcovr <- (mean(stuff$quantiles[studyp1rows, 5] > p1s & stuff$quantiles[studyp1rows, 1] < p1s))
    output[i,]$RRavgpctbias <- mean((stuff$statistics[studyRRrows, 1] - exp(lnRRs))/exp(lnRRs))
    output[i,]$p1avgpctbias <- mean((stuff$statistics[studyp1rows, 1] - p1s)/p1s)
    
    # data for the produced prior: 
    output[i,]$newRR <- stuff$statistics["RR.new", 1]
    output[i,]$newRR2point5 <- stuff$quantiles["RR.new", 1]
    output[i,]$newRR97point5 <- stuff$quantiles["RR.new", 5]
    output[i,]$newRRCIlen <- stuff$quantiles["RR.new", 5] - stuff$quantiles["RR.new", 1]
    output[i,]$newp1 <- stuff$statistics["p1.new", 1]
    output[i,]$newp1CIlen <- stuff$quantiles["p1.new", 5] - stuff$quantiles["p1.new", 1]
    output[i,]$pi11hat <- sum(n11)/sum(n)
  }
  
  avg_samps = total_samps / iter
  
  if(return_samps & ! return_counts){
    output <- list(output, avg_samps)
  }else if(!return_samps & return_counts){
    output <- list(output, list("n1" = n1, "n11" = n11))
  }else if(return_samps & return_counts){
    output <- list(output, avg_samps, list("n1" = n1, "n11" = n11))
  }
  return(output)  
}




# Test parameters ########################################################

# iter = 1; nexp = 10; n = 100; p1_alpha = 15; p1_beta = 30; lnRR_mu = 0 
# lnRR_sig = 0.1; prior_p1mu_alpha = 1; prior_p1mu_beta = 1
# prior_p1rho_alpha = .1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0
# prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2


# Plot the priors on RR and P1 for the new experiment #############################

RRsds <- c(0.05, 0.1, 0.2)
#RRsds <- c(0.3)
Ns <- c(30, 100, 300)

# overall log RR: 
lnRRmu <- log(0.7)

nexp <- 10
lnRRs <- rnorm(nexp, mean = lnRRmu, sd = 1) # Don't change this sd from 1! Change RRsds - it gets rescaled in map_comp_plot
# mu = 0.4, rho = 20
# 
p1s <- rbeta(nexp, shape1 = 8, shape2 = 12)

ess.logRR <- expand.grid(RRsds, Ns)
colnames(ess.logRR) <- c('sd', 'n')
ess.logRR$ess <- NA

t0 <- Sys.time()

for(sd in RRsds){
  for(N in Ns){
    out <- map_paper_sims(lnRRs = lnRRs, p1s = p1s, n = N, 
                         lnRRsd = sd, lnRR_mu = lnRRmu, return_samps = F, iter = 1000, 
                         return_counts = F, p1_true_mu = 0.4, p1_true_rho = 20)
    
    nexp <- length(lnRRs)
    
    out$sd <- sd
    out$N <- N
    
    if(sd == RRsds[1] & N == Ns[1]){
      gdat <- out
    }else{
      gdat <- bind_rows(gdat, out)
    }
    
    # ESS: record effective sample size for RR distributions:
    # logRRnew <- log(nextdat$RR)
    # 
    # fit <- automixfit(logRRnew, type = "norm")
    # ess.logRR[ess.logRR[,'n']==N & ess.logRR[,'sd']==sd,'ess'] <- ess(fit, sigma = 1)
  }
}

write.csv(gdat, "tblresults.csv")

print(Sys.time() - t0)
Sys.time() - t0