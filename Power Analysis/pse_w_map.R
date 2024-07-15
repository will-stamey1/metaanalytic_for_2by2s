# Probability of Successful Experiment Using MAP as Design Prior

library(rjags)
library(bayesplot)
library(ggplot2)
library(dplyr)

theme_set(theme_classic())

# Set seed 
set.seed(1234)

#  MAP JAGS model
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
  
  thetanew ~ dnorm(lnRR_mu, lnRR_tau)
  log(RRnew) <- thetanew

  theta2 ~ dbeta(p1_alpha, p1_beta)
  invRRnew <- 1/RRnew
  upboundnew <- min(1, invRRnew)
  p1new <- theta2 * upboundnew

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

", file = "map_hier")

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
  jm <- jags.model("map_hier", data = data)
  
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

# get the simulated distributions of RR* and p1*:
RRstar <- samps[[1]][, "RRnew"]
p1star <- samps[[1]][, "p1new"]

# Probability of Success JAGS model
cat(" 
model {
 
  # Prior
  
  p1 ~ dbeta(1, 1)
  RR ~ dgamma(.1, .1)
 
  p2 <- p1 * RR
    
  # likelihood
  n1 ~ dbin(p1, n)
  n11 ~ dbin(p2, n1)
  
  test <- RR - 1
  
  p.value <- 1 - step(test)

 
}

", file = "one_samp_noninf")


# Get size of generated empirical distributions:
B <- length(RRstar)

# Size of new experiment: 
n <- 500

# store results
results <- vector("numeric", B)

# iterate through logRRs and P1s from MAP simulation: 
for(i in 1:B){
  
  # Design prior:
  
  # get ith rr and ith p1:
  p1 <- p1star[i]
  rr <- RRstar[i]

  p2 <- p1 * rr
  n1 <- rbinom(1, n, p1)
  n11 <- rbinom(1, n1, p2)
  print(paste0("p1: ", p1))
  print(paste0("n1: ", n1))
  print(paste0("p2: ", p2))
  print(paste0("n11: ", n11))
  print(paste0("RR: ", rr))
  data <- list("n" = n, "n1" = n1, "n11" = n11)
  
  inits <- list(RR = 1, p1 = .5)
  
  jm <- jags.model("one_samp_noninf", data = data, inits = inits)
  samps <- update(jm, 2000)
  
  samps <- coda.samples(jm, variable.names = c("p.value"),  n.iter = 10000)
  
  stuff <- summary(samps)
  
  results[i] <- (stuff$statistics[1]>.95)
  
  rm(stuff)
  rm(samps)
  rm(jm)
}


mean(results)



