# MAP just the model 

# Load required libraries
library(rjags)
library(bayesplot)
library(ggplot2)
library(dplyr)

theme_set(theme_classic())

# Set seed 
set.seed(1234)

#  JAGS model
cat(" 
model {
 
  # Prior
  mu_p ~ dbeta(1, 1)
  rho_p ~ dgamma(.1, .1)
    
  p1_alpha <- mu_p*rho_p
  p1_beta <- rho_p * (1 - mu_p)
    
  lnRR_mu ~ dnorm(0, 0.01)
  lnRR_tau ~ dunif(0, 2)
 
  for(i in 1:nexp){  
    p2[i] <- p1[i] * RR[i] 
  }
  p2new <- p1new * RRnew
    
  # New Study Data #########
  n1new ~ dbin(p1new, n[nexp + 1])
  n11new ~ dbin(p2new, n1new)
  
  p1new ~ dbeta(p1_alpha, p1_beta)

  thetanew ~ dnorm(lnRR_mu, lnRR_tau)
  log(RRnew) <- thetanew

  # Historic Data #########
  for(i in 1:nexp){
    n1[i] ~ dbin(p1[i], n[i])
    n11[i] ~ dbin(p2[i], n1[i])
    
    p1[i] ~ dbeta(p1_alpha, p1_beta)
    
    log(RR[i]) <- theta[i]
    theta[i] ~ dnorm(lnRR_mu, lnRR_tau)
  }
}

", file = "one_samp_zero")

# PARAMETERS AND DATA INIT. ##########################

# RR parameters: 
p1_alpha <- 15
p1_beta <- 30

lnRR_mu <- -.33
lnRR_sig2 <- 0.01

nexp <- 10 # Set number of historic experiments:

n <- rep(100, nexp + 1) # Set experiment size (assume for now all experiments have same size):

iter <- 5 # simulation iterations

# output sets aside space for following values:
# 1. RRmean: mean of RR means across experiments.
# 2. RRvar: var of RR means across experiments. 
# 3. RRCIlen: Average CI length for RR.
# 3. RRmeancovr: coverage rate of RR central mean 
# 4. p1meancovr: coverage rate of p1 central mean
# 5. RRindivcovr: coverage rate for experiment RR means
# 6. p1indivcovr: coverage rate for experiment p1 means
# 7. newRR: mean of new experiment's RR mcmc chain
# 8. newRRCIlen: CI length from new experiment's RR mcmc chain.

output <- data.frame(matrix(NA, iter, 11))
colnames(output) <- c("RRbar", "RRvar", "RRCIlen", "RRmeancovr",
                      "p1meancovr", "RRindivcovr", "p1indivcovr", 
                      "newRR", "newRRCIlen", "newp1", "newp1CIlen")

for(i in 1:iter){
  
  # sample a p1 and RR for each historic experiment: 
  p1 <- rbeta(n = nexp, p1_alpha, p1_beta)  
  lnRR <- rnorm(nexp, mean = lnRR_mu, lnRR_sig2)
  p1_trumu <- p1_alpha/(p1_alpha + p1_beta)
  
  # Generate data according to the parameters: 
  n1 <- rbinom(nexp, n, prob = p1)
  n11 <- rbinom(nexp, n1, prob = (p1 * exp(lnRR)))
  
  # Define data for JAGS
  data <- list(n = n, nexp = nexp, n1 = n1, n11 = n11)
  
  # Compile
  jm <- jags.model("one_samp_zero", data = data)
  
  # Burnin
  samps <- update(jm, 5000)
  
  # iterations
  samps <- coda.samples(jm, variable.names = c("p1", "RR", "p1new", "RRnew"),  n.iter = 20000)
  
  stuff <- summary(samps)
  
  output[i,]$RRbar <- mean(stuff$statistics[1:nexp,1])
  output[i,]$RRvar <- var(stuff$statistics[1:nexp,1])
  output[i,]$RRCIlen <- mean(stuff$quantiles[1:nexp,5]) - mean(stuff$quantiles[,1])
  output[i,]$RRmeancovr <- (mean(stuff$quantiles[1:nexp,5]) > exp(lnRR_mu) & mean(stuff$quantiles[1:nexp,1]) < exp(lnRR_mu))
  output[i,]$p1meancovr <- (mean(stuff$quantiles[nexp + 1:nexp, 5]) > p1_trumu & mean(stuff$quantiles[nexp + 1:nexp, 1]) < p1_trumu)
  output[i,]$RRindivcovr <- (mean(stuff$quantiles[1:nexp,5] > exp(lnRR) & stuff$quantiles[1:nexp,1] < exp(lnRR)))
  output[i,]$p1indivcovr <- (mean(stuff$quantiles[nexp + 1:nexp, 5] > p1 & stuff$quantiles[nexp + 1:nexp, 1] < p1))
  output[i,]$newRR <- stuff$statistics["RRnew", 1]
  output[i,]$newRRCIlen <- stuff$quantiles["RRnew", 5] - stuff$quantiles["RRnew", 1]
  output[i,]$newp1 <- stuff$statistics["p1new", 1]
  output[i,]$newp1CIlen <- stuff$quantiles["p1new", 5] - stuff$quantiles["p1new", 1]
}

mean(output1[,1])
mean(output1[,2])
mean(output1[,3])

# Plot the priors on RR and P1 for the new experiment #############################

gdat <- data.frame(samps[[1]][,"RRnew"],
                   samps[[1]][,"p1new"])
colnames(gdat) <- c("RR", "p1")

# RR
bayesplot::mcmc_trace(samps, pars = "RRnew")
bayesplot::mcmc_hist(samps, pars = "RRnew", binwidth = 0.5)

gdat %>% ggplot(aes(x = RR)) + 
  geom_histogram(binwidth = 0.05) + 
  scale_x_continuous(limits = c(0, 6))

# P1
samps[[1]][,"p1new"]
bayesplot::mcmc_trace(samps, pars = "p1new")
bayesplot::mcmc_hist(samps, pars = "p1new", binwidth = 0.01)

gdat %>% ggplot(aes(x = p1)) + 
  geom_histogram(binwidth = 0.01) + 
  scale_x_continuous(limits = c(0, 1))
