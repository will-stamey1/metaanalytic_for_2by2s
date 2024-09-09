# Load required libraries
library(rjags)

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

lnRR_mu <- 0
lnRR_sig2 <- 0.5

nexp <- 10 # Set number of historic experiments:

n <- rep(100, nexp) # Set experiment size (assume for now all experiments have same size):

iter <- 5 # simulation iterations

# output sets aside space for following values:
# 1. expRRbar: mean of exp RR means across experiments.
# 2. expRRvar: var of exp RR means across experiments. 
# 3. RRCIlen: Average CI length for RR.
# 3. RRmeancovr: coverage rate of RR central mean 
# 4. p1meancovr: coverage rate of p1 central mean
# 5. RRindivcovr: coverage rate for experiment RR means
# 6. p1indivcovr: coverage rate for experiment p1 means
output <- data.frame(matrix(NA, iter, 7))
colnames(output) <- c("expRRbar", "expRRvar", "RRCIlen", "RRmeancovr",
                      "p1meancovr", "RRindivcovr", "p1indivcovr")

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
  samps <- coda.samples(jm, variable.names = c("p1", "RR", "mu_p", "rho_p"),  n.iter = 20000)
  
  stuff <- summary(samps)
  
  output[i,]$expRRbar <- mean(stuff$statistics[1:nexp,1])
  output[i,]$expRRvar <- var(stuff$statistics[1:nexp,1])
  output[i,]$RRCIlen <- mean(stuff$quantiles[1:nexp,5]) - mean(stuff$quantiles[,1])
  output[i,]$RRmeancovr <- (mean(stuff$quantiles[1:nexp,5]) > exp(lnRR_mu) & mean(stuff$quantiles[1:nexp,1]) < exp(lnRR_mu))
  output[i,]$p1meancovr <- (mean(stuff$quantiles[nexp + 1:nexp, 5]) > p1_trumu & mean(stuff$quantiles[nexp + 1:nexp, 1]) < p1_trumu)
  output[i,]$RRindivcovr <- (mean(stuff$quantiles[1:nexp,5] > exp(lnRR) & stuff$quantiles[1:nexp,1] < exp(lnRR)))
  output[i,]$p1indivcovr <- (mean(stuff$quantiles[nexp + 1:nexp, 5] > p1 & stuff$quantiles[nexp + 1:nexp, 1] < p1))
}

mean(output1[,1])
mean(output1[,2])
mean(output1[,3])


