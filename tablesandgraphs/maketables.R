# Analyze tablesandgraphs simulation results:

library(dplyr)
library(stargazer)
library(tidyr)

a <- read.csv("tablesandgraphs/table1and2job/tblresults8-14-24_unifsdprior.csv")
#a <- read.csv("tablesandgraphs/table1and2job/tblresults8-09-24_.55.csv")

# FORGOT TO ADD LABELS FOR lnRRmu, ADD MANUALLY:

# a$lnRRmu <- NA
# a[1001:8000,]$lnRRmu <- 0.1823216
# a[15001:16000,]$lnRRmu <- 0.1823216
# a[is.na(a["lnRRmu"]),]$lnRRmu <- -0.3566749

# Calculate avg pct bias for parameters: 
a <- 
  a %>% 
  mutate(
    lnRRmu_pctbias = (lnRRmu_mean - lnRRmu)/lnRRmu * 100,
    lnRRsd_pctbias = (lnRRsd_mean - sd)/sd * 100,
    mup_pctbias = (mup_mean - 0.4)/0.4 * 100, # This of course must change if mup is changed in the sims
    rhop_pctbias = (rhop_mean - 20)/20 * 100 # This must change if rhop is changed in the sims
  ) %>% 
  relocate(lnRRmu_pctbias, .after=lnRRmu_mean) %>% 
  relocate(lnRRsd_pctbias, .after=lnRRsd_mean) %>% 
  relocate(mup_pctbias, .after=mup_mean) %>% 
  relocate(rhop_pctbias, .after=rhop_mean) %>% 
  select(!contains("_mean")) %>% # drop mean variables 
  # multiply RRi and pi pct bias by a hundred, calculated this wrong originally: 
  mutate(RRavgpctbias  = RRavgpctbias * 100,
         p1avgpctbias = p1avgpctbias * 100) %>% 
  relocate(RRavgpctbias, .before=RRindivcovr) %>% 
  relocate(p1avgpctbias, .before=p1indivcovr)
  


# Table 1: 

t1 <- 
  a %>% group_by(nexp, sd, N, exp(lnRRmu)) %>% 
  summarize(
    #mean(exp(lnRRmu_mean)),
    across(lnRRmu_pctbias:p1indivcovr, ~mean(.x))) %>% 
  arrange(sd, nexp, N, `exp(lnRRmu)`) %>% 
  relocate(sd, nexp, N)

# Table 2:

t2 <- a %>% 
  group_by(nexp, sd, N, exp(lnRRmu)) %>% 
  summarize(
    across(newRR:newp1plusCIlen, ~mean(.x))
  ) %>% 
  arrange(sd, nexp, N, `exp(lnRRmu)`) 



# Make latex output ###########################################################

# table 1
t1 %>%
  #select(contains("exp(lnRRmu)") | (!contains("lnRRmu") & !contains("lnRRsd"))) %>% 
  select(c("sd", "nexp", "N", "exp(lnRRmu)", "lnRRmu_pctbias", "lnRRmu_CIlen", "lnRRmu_covr", "lnRRsd_pctbias", "lnRRsd_CIlen", "lnRRsd_covr")) %>% 
  relocate(sd, nexp, N) %>% 
  mutate(across(everything(), ~round(.x,3))) %>% 
  #pivot_longer(cols = lnRRmu_mean:p1avgpctbias) %>% 
  arrange(`exp(lnRRmu)`) %>% 
  stargazer(type = "latex", header = T, summary = F, vars = colnames(t1), colnames = T)

# table 2
t1 %>%
  select(contains("exp(lnRRmu)") | (!contains("lnRRmu") & !contains("lnRRsd"))) %>% 
  select(1:10) %>% 
  relocate(sd, nexp, N) %>% 
  mutate(across(everything(), ~round(.x,3))) %>% 
  #pivot_longer(cols = lnRRmu_mean:p1avgpctbias) %>% 
  arrange(`exp(lnRRmu)`) %>% 
  stargazer(type = "latex", header = T, summary = F, vars = colnames(t1), colnames = T)

# table 3
t1 %>%
  select(contains("exp(lnRRmu)") | (!contains("lnRRmu") & !contains("lnRRsd"))) %>% 
  select(c(1:4, 11:14)) %>% 
  relocate(sd, nexp, N) %>% 
  mutate(across(everything(), ~round(.x,3))) %>% 
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
