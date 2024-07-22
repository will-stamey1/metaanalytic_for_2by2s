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
  
  q1 ~ dbeta(30.88, 49.12)
  # q1 ~ dbeta(1, 1)
  RR2 <- 1/RR
  p1 <- min(q1, RR2)
  RR ~ dgamma(.1, .1)
 
  
    
  # likelihood
  z[1:3] ~ dmulti(p[],n)
  p[1]<- RR*pow(p1,2)
	p[2]<- p1 - RR*pow(p1,2)
	p[3]<- 1 - p1
  
  test <- RR - 1
  
  p.value <- 1 - step(test)

 
}

", file = "one_samp_noninf")






#simulation size
B <- 1000

#sample size
n <- 280

# store results
output <- vector("numeric", B)

for(i in 1:B) {

  # Design prior
  
  logrr <- rnorm(1, log(.7), .1)
  rr <- .576 #exp(logrr)
  p1 <- .386 #rbeta(1, 38.6, 61.4)
  p1 <- min(p1, 1/rr)
  
  pa <- rr*p1^2
  pb <- p1 - rr*p1^2
  pc <- 1 - p1
  
  p <- cbind(pa, pb, pc)
  
  z <- rmultinom(1, n, p)
  
  data <- list("n" = n, z = as.numeric(z))
  
  inits <- list(RR = 1, q1 = .5)
  
  jm <- jags.model("one_samp_noninf", data = data, inits = inits)
  samps <- update(jm, 2000)
  
  samps <- coda.samples(jm, variable.names = c( "p.value"),  n.iter = 10000)
  
  stuff <- summary(samps)
  
  output[i] <- (stuff$statistics[1] > .975)
}


mean(output)



