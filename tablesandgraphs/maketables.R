# Analyze tablesandgraphs simulation results:

library(dplyr)

a <- read.csv("tablesandgraphs/table1and2job/tblresults.csv")

a$lnRRmu_mean %>% exp() %>%  hist()

mean(exp(a$lnRRmu_mean))

a$mup_mean %>% mean()
a$lnRRsd_covr %>% mean()

