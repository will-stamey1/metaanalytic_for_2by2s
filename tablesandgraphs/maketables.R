# Analyze tablesandgraphs simulation results:

library(dplyr)
library(stargazer)
library(tidyr)

a <- read.csv("tablesandgraphs/table1and2job/tblresults8-09-24.csv")

a$lnRRmu_mean %>% exp() %>% hist()

mean(exp(a$lnRRmu_mean))

a$mup_mean %>% mean()
a$lnRRsd_covr %>% mean()

# FORGOT TO ADD LABELS FOR lnRRmu, ADD MANUALLY:

a$lnRRmu <- NA
a[1001:8000,]$lnRRmu <- 0.1823216
a[15001:16000,]$lnRRmu <- 0.1823216
a[is.na(a["lnRRmu"]),]$lnRRmu <- -0.3566749

# Table 1: 

t1 <- 
  a %>% group_by(nexp, sd, N, exp(lnRRmu)) %>% 
  summarize(
    #mean(exp(lnRRmu_mean)),
    across(lnRRmu_mean:p1avgpctbias, ~mean(.x))) %>% 
  arrange(sd, nexp, N, `exp(lnRRmu)`) %>% 
  relocate(nexp, N, sd)

# Table 2:

t2 <- a %>% 
  group_by(nexp, sd, N, exp(lnRRmu)) %>% 
  summarize(
    across(newRR:newp1plusCIlen, ~mean(.x))
  ) %>% 
  arrange(sd, nexp, N, `exp(lnRRmu)`) %>% 
  relocate(nexp, N, sd)



# Make latex output ###########################################################

t1 %>% mutate(across(everything(), ~round(.x,3))) %>% 
  #pivot_longer(cols = lnRRmu_mean:p1avgpctbias) %>% 
  arrange(`exp(lnRRmu)`) %>% 
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
# Save tables to csvs #########################################################

write.csv(t1, "t1.csv")
write.csv(t2, "t2.csv")
