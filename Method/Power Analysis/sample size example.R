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






#simulation size
B <- 100

#sample size
n <- 500

# store results
output <- vector("numeric", B)

for(i in 1:B) {

# Design prior


logrr <- rnorm(1, log(.7), .1)
rr <- .7 #exp(logrr)
p1 <- rbeta(1, 10, 30)
p1 <- min(p1, 1/rr)

p2 <- p1 * rr
n1 <- rbinom(1, n, p1)
n11 <- rbinom(1, n1, p2)

data <- list("n" = n, "n1" = n1, "n11" = n11)

inits <- list(RR = 1, p1 = .5)

jm <- jags.model("one_samp_noninf", data = data, inits = inits)
samps <- update(jm, 2000)

samps <- coda.samples(jm, variable.names = c( "p.value"),  n.iter = 10000)

stuff <- summary(samps)

output[i] <- (stuff$statistics[1]>.95)

}


mean(output)



