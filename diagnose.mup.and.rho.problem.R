# figuring out what in the world is wrong with a simple demo: 

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

expit <- function(x){
  exp(x)/(1 + exp(x))
}

# Test parameters ########################################################

# iter = 1; nexp = 10; n = 100; p1_alpha = 15; p1_beta = 30; lnRR_mu = 0 
# lnRR_sig = 0.1; prior_p1mu_alpha = 1; prior_p1mu_beta = 1
# prior_p1rho_alpha = .1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0
# prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2


# Set parameters
RRsds <- c(0.05, 0.1, 0.2)
#RRsds <- c(0.3)
Ns <- c(30, 100)
nexps <- c(5, 15, 25)

n = Ns[1] # number of samples per experiment
lnRRsd = RRsds[3] # true lnRR sd 
nexp <- nexps[3]
lnRR_mu = log(0.7) # true lnRR mu
p1_true_mu = 0.4
p1_true_rho = 20
iter = 1
  
# Prior parameters: 
prior_p1mu_alpha = 1; prior_p1mu_beta = 1; 
prior_p1rho_alpha = 1; prior_p1rho_beta = .1; prior_lnRRmu_mu = 0; 
prior_lnRRmu_tau = .01; prior_lnRRtau_lo = 0; prior_lnRRtau_hi = 2;

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
    p11[i]<- RR[i]*pow(p1plus[i],2)
    p12[i]<- p1plus[i] - RR[i]*pow(p1plus[i],2)
    p22[i]<- 1 - p1plus[i]
    
    ps[i,1] <- p11[i]
    ps[i,2] <- p12[i]
    ps[i,3] <- p22[i]
    z[i,1:3] ~ dmulti(ps[i,1:3],n[i])

    theta1[i] ~ dbeta(p1_alpha, p1_beta)
    invRR[i] <- 1/RR[i]
    upbound[i] <- min(1, invRR[i]) 
    p1plus[i] <- theta1[i] * upbound[i] 
    
    log(RR[i]) <- theta[i]
    theta[i] ~ dnorm(lnRR_mu, lnRR_tau)
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
output <- data.frame(matrix(1, iter, 21))
colnames(output) <- c("lnRRmu_mean", "lnRRmu_CIlen", "lnRRmu_covr", 
                      "lnRRsd_mean", "lnRRsd_CIlen", "lnRRsd_covr",
                      "mup_mean", "mup_CIlen", "mup_covr", 
                      "rhop_mean", "rhop_CIlen", "rhop_covr",
                      "RRindivcovr", "RRavgpctbias", "p1indivcovr", "p1avgpctbias",
                      "newRR", "newRRCIlen", 
                      "newp1plus", "newp1plusCIlen", "p11hat")

#for(i in 1:iter){
  
# Generate data: 
lnRRs <- rnorm(nexp, mean = lnRR_mu, sd = lnRRsd)
rrs <- exp(lnRRs)

# mu = 0.4, rho = 20
p1s <- rbeta(nexp, shape1 = 8, shape2 = 12)
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

i <- 1

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



