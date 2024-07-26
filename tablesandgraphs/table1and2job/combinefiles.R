library(dplyr)

matching_files <- list.files(path = "", pattern = "tbl1and2sims_")

dfs <- lapply(matching_files, read.csv)

fulldf <- bind_rows(dfs)

write.csv(fulldf, "tblresults7-26-24.csv")