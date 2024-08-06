# Analyze tablesandgraphs simulation results:

library(dplyr)
library(stargazer)

a <- read.csv("tablesandgraphs/table1and2job/tblresults7-26-24_2.csv")

a$lnRRmu_mean %>% exp() %>%  hist()

mean(exp(a$lnRRmu_mean))

a$mup_mean %>% mean()
a$lnRRsd_covr %>% mean()

# FORGOT TO ADD LABELS FOR CONFIGURATION PARAMETERS, ADD MANUALLY:
# a$sd <- NA
# a$N <- NA
# 
# RRsds <- c(0.05, 0.1, 0.2)
# Ns <- c(30, 100, 300)
# 
# i <- 0
# for(sd in RRsds){
#   for(N in Ns){
#     
#     a[(1000*i + 1):(1000*(i + 1)),]$sd <- sd
#     a[(1000*i + 1):(1000*(i + 1)),]$N <- N
#     i <- i + 1
#   }
# }

# Table 1: 

t1 <- 
  a %>% group_by(sd, N) %>% 
  summarize(
    #mean(exp(lnRRmu_mean)),
    across(lnRRmu_mean:p1avgpctbias, ~mean(.x)))

# Table 2:

t2 <- a %>% 
  group_by(sd, N) %>% 
  summarize(
    across(newRR:newp1plusCIlen, ~mean(.x))
  )


# Make latex output ###########################################################

t1 %>% mutate(across(everything(), ~round(.x,3))) %>% 
  stargazer(type = "latex", header = T, summary = F, vars = colnames(t1), colnames = T)

t2 %>% mutate(across(everything(), ~round(.x,3))) %>% 
  stargazer(type = "latex", header = T, summary = F, vars = colnames(t1), colnames = T)



# Histograms ####################################################################

asd01 <- a %>% filter(sd == 0.1, N == 100) 
hist(asd01$rhop_mean)

a %>% filter(sd == 0.1) %>% 
  mutate(N = factor(N)) %>% 
  ggplot(aes(x = rhop_mean)) +
  geom_histogram(group = N, fill = N)
