# Load required libraries
library(rjags)

# Set seed 
set.seed(1234)

# parameters
n <- 100  # Number of trials
p1 <- .3
RR <- 1.2



#  JAGS model
cat(" 
model {

 n1 ~ dbin(p1, n)
 
 n11 ~ dbin(p2, n1)
 
  p2 <- p1*RR
 
 # p2 ~ dbeta(1, 1)

  # Prior
  p1 ~ dbeta(1, 1)

  log(RR) <- theta
 
 theta ~ dnorm(0, .1)
 
# RR <- p2/p1
  
}

", file = "one_samp_zero")

B <- 500
n <- 300


output1 <- matrix(NA, B, 3)

for(i in 1:B) {
  

  n1 <- rbinom(1, n, p1)  
  n11 <- rbinom(1, n1, p1*RR)
  
  # Define data for JAGS
  data <- list(n = n, n1 = n1, n11 = n11)
  
  
  # Compile
  jm <- jags.model("one_samp_zero", data = data)
  
  # Burnin
  samps <- update(jm, 5000)
  
  # iterations
  samps <- coda.samples(jm, variable.names = c( "RR"),  n.iter = 20000)
  
  
  stuff <- summary(samps)
  
  output1[i, 1] <- stuff$statistics[1]
  output1[i, 2] <- stuff$quantiles[5] - stuff$quantiles[1]
  output1[i, 3] <- (stuff$quantiles[5] > RR & stuff$quantiles[1] < RR)

}

mean(output1[,1])
mean(output1[,2])
mean(output1[,3])


