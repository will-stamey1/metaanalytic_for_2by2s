# Load required libraries
library(rjags)

# Set seed 
set.seed(1234)

hier2b2_logistp <- function(iter = 1, nexp = 10, n = 100, p1_alpha = 15, p1_beta = 30, lnRR_mu = 0, 
                    lnRR_sig = 0.1, prior_p1mu_mu = 0, prior_p1mu_tau = .01,
                    prior_p1_sig_lo = 0, prior_p1_sig_hi = 2, prior_lnRRmu_mu = 0, 
                    prior_lnRRmu_tau = .01, prior_lnRRtau_lo = 0, prior_lnRRtau_hi = 2){
  
  # nexp: number of experiments in the simulation.
  # iter: simulation iterations.
  # p1_alpha: True alpha value of the population p1.
  # p1_beta: True beta value of the population p1.
  # lnRR_mu: True mu value of population lnRR.
  # lnRR_sig: True std. dev. of population lnRR.
  # prior_p1mu_mu and prior_p1mu_tau: parameters for normal prior on p1s mean, p1mu.
  # prior_p1_tau_lo and prior_p1_tau_hi: Lower and upper bounds for the uniform prior of the variance of logit(p)
  # prior_lnRRmu_mu and prior_lnRRmu_tau: parameters for normal prior of lnRRmu.
  # prior_lnRRtau_lo and prior_lnRRtau_hi: parameters for uniform prior for lnRRtau.
  
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
    p1_mu ~ dnorm(prior_p1mu_mu, prior_p1mu_tau)
    p1_sig ~ dunif(prior_p1_sig_lo, prior_p1_sig_hi)
    p1_tau <- 1/(p1_sig*p1_sig)  
      
    lnRR_mu ~ dnorm(prior_lnRRmu_mu, prior_lnRRmu_tau)
    lnRR_sig ~ dunif(prior_lnRRtau_lo, prior_lnRRtau_hi)
    lnRR_tau <- 1/(lnRR_sig*lnRR_sig)
   
    for(i in 1:nexp){  
      p2[i] <- p1[i] * RR[i] 
    }
      
    for(i in 1:nexp){
      n1[i] ~ dbin(p1[i], n[i])
      n11[i] ~ dbin(p2[i], n1[i])
      
      logitp1[i] ~ dnorm(p1_mu, p1_tau)
      logit(p1[i]) <- logitp1[i]
      
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
  # 4. lnRRtau_mean: mean of lnRR_tau chain. 
  # 5. lnRRtau_CIlen: CI length for lnRR_tau.
  # 6. lnRRtau_covr: coverage of lnRR_tau true value.
  # 7. mup_mean: mean of p1 mu chain.
  # 8. mup_CIlen: CI length for p1 mu chain. 
  # 9. mup_covr: coverage rate of p1 central mean
  # 10. rhop_mean: mean of p1 rho.
  # 11. rhop_CIlen: CI length for p1 rho chain.
  # 12. rhop_covr: coverage of p1 rho true value. . 
  # 13. RRindivcovr: coverage rate for experiment RR means
  # 14. p1indivcovr: coverage rate for experiment p1 means
  output <- data.frame(matrix(NA, iter, 15))
  colnames(output) <- c("lnRRmu_mean", "lnRRmu_CIlen", "lnRRmu_covr", 
                        "lnRRsig_mean", "lnRRsig_CIlen", "lnRRsig_covr",
                        "p1mu_mean", "p1mu_CIlen", "p1mu_covr", 
                        "p1sig_mean", "p1sig_CIlen", "p1sig_covr",
                        "RRindivcovr", "p1indivcovr", "p1pctbias")
  
  for(i in 1:iter){
    
    # sample a p1 and RR for each historic experiment: 
    p1 <- rbeta(n = nexp, p1_alpha, p1_beta)  
    lnRR <- rnorm(nexp, mean = lnRR_mu, lnRR_sig)
    p1_true_mu <- p1_alpha/(p1_alpha + p1_beta)
    p1_true_rho <- p1_alpha + p1_beta
    
    # Generate data according to the parameters: 
    n1 <- rbinom(nexp, n, prob = p1)
    n11 <- rbinom(nexp, n1, prob = (p1 * exp(lnRR)))
    
    # Define data for JAGS
    data <- list(n = n, nexp = nexp, n1 = n1, n11 = n11, 
                 prior_p1mu_mu = prior_p1mu_mu, 
                 prior_p1mu_tau = prior_p1mu_tau, prior_p1_sig_lo = prior_p1_sig_lo, 
                 prior_p1_sig_hi = prior_p1_sig_hi,
                 prior_lnRRmu_mu = prior_lnRRmu_mu, prior_lnRRmu_tau = prior_lnRRmu_tau,
                 prior_lnRRtau_lo = prior_lnRRtau_lo, prior_lnRRtau_hi = prior_lnRRtau_hi)
    
    # Compile
    jm <- jags.model("one_samp_zero", data = data)
    
    # Burnin
    samps <- update(jm, 5000)
    
    
    # iterations
    samps <- coda.samples(jm, variable.names = c("lnRR_mu", "lnRR_sig", "p1_mu", 
                                                 "p1_sig", "p1", "RR"),  
                          n.iter = 20000)
    
    stuff <- summary(samps)
    
    output[i,]$lnRRmu_mean <- stuff$statistics["lnRR_mu",1]
    output[i,]$lnRRmu_CIlen <- stuff$quantiles["lnRR_mu",5] - stuff$quantiles["lnRR_mu",1]
    output[i,]$lnRRmu_covr <- stuff$quantiles["lnRR_mu",5] > lnRR_mu & stuff$quantiles["lnRR_mu",1] < lnRR_mu
    output[i,]$lnRRsig_mean <- stuff$statistics["lnRR_sig",1]
    output[i,]$lnRRsig_CIlen <- stuff$quantiles["lnRR_sig",5] - stuff$quantiles["lnRR_sig",1]
    output[i,]$lnRRsig_covr <- stuff$quantiles["lnRR_sig",5] > sqrt(lnRR_sig) & stuff$quantiles["lnRR_sig",1] < sqrt(lnRR_sig)
    output[i,]$p1mu_mean <- stuff$statistics["p1_mu",1]
    #output[i,]$p1mu_CIlen <- stuff$quantiles["p1_mu",5] - stuff$quantiles["p1_mu",1]
    #output[i,]$p1mu_covr <- stuff$quantiles["p1_mu",5] > p1_true_mu & stuff$quantiles["p1_mu",1] < p1_true_mu
    output[i,]$p1sig_mean <- stuff$statistics["p1_sig",1]
    #output[i,]$p1sig_CIlen <- stuff$quantiles["p1_sig",5] - stuff$quantiles["p1_sig",1]
    #output[i,]$p1sig_covr <- stuff$quantiles["p1_sig",5] > p1_true_rho & stuff$quantiles["p1_sig",1] < p1_true_rho
    output[i,]$RRindivcovr <- (mean(stuff$quantiles[1:nexp,5] > exp(lnRR) & stuff$quantiles[1:nexp,1] < exp(lnRR)))
    output[i,]$p1indivcovr <- (mean(stuff$quantiles[(nexp + 2) + 1:nexp, 5] > p1 & stuff$quantiles[(nexp + 2) + 1:nexp, 1] < p1))
    output[i,]$p1pctbias <- mean((p1 - stuff$statistics[grep("p1\\[", rownames(stuff$statistics))])/p1)
  }
  
  return(output)
  
}


hier2b2_logistp(iter = 10)

