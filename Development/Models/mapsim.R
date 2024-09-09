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

mapsim <- function(iter = 1, n.new = n, RR.new = -.4, p1.new = .3,
                   nexp = 10, n = 30, p1_alpha = 15, p1_beta = 30, lnRR_mu = -.3, 
                    lnRR_sig = 0.1, prior_p1mu_alpha = 1, prior_p1mu_beta = 1, 
                    prior_p1rho_alpha = .1, prior_p1rho_beta = .1, prior_lnRRmu_mu = 0, 
                    prior_lnRRmu_tau = .01, prior_lnRRtau_lo = 0, prior_lnRRtau_hi = 2,
                    return_samps = F){

  # nexp: number of experiments in the simulation.
  # iter: simulation iterations.
  # p1_alpha: True alpha value of the population p1.
  # p1_beta: True beta value of the population p1. 
  # lnRR_mu: True mu value of population lnRR. 
  # lnRR_sig2: True variance of population lnRR.
  # prior_p1mu_alpha and prior_p1mu_beta: parameters for beta prior on p1
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
  output <- data.frame(matrix(NA, iter, 19))
  colnames(output) <- c("lnRRmu_mean", "lnRRmu_CIlen", "lnRRmu_covr", 
                        "lnRRtau_mean", "lnRRtau_CIlen", "lnRRtau_covr",
                        "mup_mean", "mup_CIlen", "mup_covr", 
                        "rhop_mean", "rhop_CIlen", "rhop_covr",
                        "RRindivcovr", "p1indivcovr",
                        "newRR", "newRRCIlen", "newp1", "newp1CIlen",
                        "pi11hat")
  
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
    samps <- coda.samples(jm, variable.names = c("lnRR_mu", "lnRR_tau", "mu_p", 
                                                 "rho_p", "p1", "RR",
                                                 "p1.new", "RR.new"),  
                          n.iter = 20000)
    
    stuff <- summary(samps)
    
    output[i,]$lnRRmu_mean <- stuff$statistics["lnRR_mu",1]
    output[i,]$lnRRmu_CIlen <- stuff$quantiles["lnRR_mu",5] - stuff$quantiles["lnRR_mu",1]
    output[i,]$lnRRmu_covr <- stuff$quantiles["lnRR_mu",5] > lnRR_mu & stuff$quantiles["lnRR_mu",1] < lnRR_mu
    output[i,]$lnRRtau_mean <- stuff$statistics["lnRR_tau",1]
    output[i,]$lnRRtau_CIlen <- stuff$quantiles["lnRR_tau",5] - stuff$quantiles["lnRR_tau",1]
    output[i,]$lnRRtau_covr <- stuff$quantiles["lnRR_tau",5] > exp(lnRR_mu) & stuff$quantiles["lnRR_tau",1] < exp(lnRR_mu)
    output[i,]$mup_mean <- stuff$statistics["mu_p",1]
    output[i,]$mup_CIlen <- stuff$quantiles["mu_p",5] - stuff$quantiles["mu_p",1]
    output[i,]$mup_covr <- stuff$quantiles["mu_p",5] > p1_true_mu & stuff$quantiles["mu_p",1] < p1_true_mu
    output[i,]$rhop_mean <- stuff$statistics["rho_p",1]
    output[i,]$rhop_CIlen <- stuff$quantiles["rho_p",5] - stuff$quantiles["rho_p",1]
    output[i,]$rhop_covr <- stuff$quantiles["rho_p",5] > p1_true_rho & stuff$quantiles["rho_p",1] < p1_true_rho
    output[i,]$RRindivcovr <- (mean(stuff$quantiles[1:nexp,5] > exp(lnRR) & stuff$quantiles[1:nexp,1] < exp(lnRR)))
    output[i,]$p1indivcovr <- (mean(stuff$quantiles[nexp + 1:nexp, 5] > p1 & stuff$quantiles[nexp + 1:nexp, 1] < p1))
    output[i,]$newRR <- stuff$statistics["RR.new", 1]
    output[i,]$newRRCIlen <- stuff$quantiles["RR.new", 5] - stuff$quantiles["RR.new", 1]
    output[i,]$newp1 <- stuff$statistics["p1.new", 1]
    output[i,]$newp1CIlen <- stuff$quantiles["p1.new", 5] - stuff$quantiles["p1.new", 1]
    output[i,]$pi11hat <- sum(n11)/sum(n)
  }
  
  if(return_samps){
    output <- list(output, samps)
  }
  return(output)  
}


a <- exp(rnorm(10000000, mean = lnRR_mu, lnRR_sig))

outpoint1 <- mapsim(iter = 5, lnRR_sig=0.1)
outpoint5 <- mapsim(iter = 5, lnRR_sig=0.5)


# ESS #######################################################

a<-mapsim(iter = 1, return_samps = T, n = 300, lnRR_sig = 0.3)

RRnew <- a[[2]][, "RR.new"]
logRRnew <- log(RRnew[[1]][1:length(RRnew[[1]])])

pi11hat <- a[[1]]$pi11hat
sigma <- sqrt((1 - pi11hat)/pi11hat)

fit <- automixfit(logRRnew, type = "norm")

ess(fit, sigma = sigma)





# Test parameters ########################################################

# iter = 1; nexp = 10; n = 100; p1_alpha = 15; p1_beta = 30; lnRR_mu = 0 
# lnRR_sig = 0.1; prior_p1mu_alpha = 1; prior_p1mu_beta = 1
# prior_p1rho_alpha = .1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0
# prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2


# Plot the priors on RR and P1 for the new experiment #############################

RRsds <- c(0.1, 0.2, 0.3)
Ns <- c(30, 100, 300)


ess.logRR <- expand.grid(RRsds, Ns)
colnames(ess.logRR) <- c('sd', 'n')
ess.logRR$ess <- NA

for(sd in RRsds){
  for(N in Ns){
    out <- mapsim(n = N, lnRR_sig = sd, return_samps = T)
    
    
    samps <- out[[2]]
    
    nextdat <- data.frame(samps[[1]][,"RR.new"],
                          samps[[1]][,"p1.new"])
    
    colnames(nextdat) <- c("RR", "p1")
    
    nextdat$sd <- sd
    nextdat$n <- N
    
    if(sd == RRsds[1] & N == Ns[1]){
      gdat <- nextdat
    }else{
      gdat <- bind_rows(gdat, nextdat)
    }
    
    # ESS: record effective sample size for RR distributions:
    # logRRnew <- log(nextdat$RR)
    # 
    # fit <- automixfit(logRRnew, type = "norm")
    # ess.logRR[ess.logRR[,'n']==N & ess.logRR[,'sd']==sd,'ess'] <- ess(fit, sigma = 1)
  }
}

sdlabs <- sdnames <- c()
for(i in unique(gdat$sd)){
  sdnames <- c(sdlabs, i)
  sdlabs <- c(sdnames, "expression(sigma)")
}
names(sdlabs) <- sdnames

# RR
gdat %>% 
  mutate(sigma = as.character(sd)) %>% 
  ggplot(aes(x = RR)) + 
  geom_histogram(binwidth = 0.05, color = "black", fill = "blue") + 
  scale_x_continuous(limits = c(0, 2)) + 
  facet_grid(sigma ~ n, labeller = labeller(.rows = label_both, .cols = label_both)) + 
  ylab("") + 
  xlab("Risk-Ratio")

ggsave("RRpriors.jpg", width = 6, height = 4)




p1_ps <- list(c(5, 10), c(10, 20), c(100, 200))
Ns <- c(30, 100, 300)

for(p1p in p1_ps){
  #p1p <- p1p[[1]]
  for(N in Ns){
    p1a <- p1p[1]
    p1b <- p1p[2]
    out <- mapsim(n = N, p1_alpha = p1a, p1_beta = p1b, return_samps = T)
    samps <- out[[2]]
    
    nextdat <- data.frame(samps[[1]][,"RR.new"],
                          samps[[1]][,"p1.new"])
    
    colnames(nextdat) <- c("RR", "p1")
    
    nextdat$p1p <- paste0(as.character(p1a), ", ", as.character(p1b))
    nextdat$n <- N
    
    if(p1a == p1_ps[[1]][1] & N == Ns[1]){
      gdat <- nextdat
    }else{
      gdat <- bind_rows(gdat, nextdat)
    }
  }
}

# make categories for beta params:
p1pnames <- c()
for(i in p1_ps){
  #j <- i[[1]]
  p1pnames <- c(p1pnames, paste0(as.character(i[1]), ", ", as.character(i[2])))
}

gdat$p1p <- factor(gdat$p1p, levels = p1pnames)

betafy <- function(x){
  # takes a string with two numbers and adds 
  # labels "alpha" and "beta"
  ab <- x %>% str_split(",")
  a <- ab[[1]][1]
  b <- ab[[1]][2]
  
  paste0("alpha: ", a, ", ", "beta: ", b)
}

# P1
gdat %>% 
  ggplot(aes(x = p1)) + 
  geom_histogram(binwidth = 0.03, color = "black", fill = "blue") + 
  scale_x_continuous(limits = c(0, 1)) + 
  facet_grid(p1p ~ n) + 
  labs(y = "", x = expression(bold(p["1"])))

