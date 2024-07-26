library(dplyr)

matching_files <- list.files(pattern = "tbl1and2sims_")

print(matching_files)

dfs <- lapply(matching_files, read.csv)

fulldf <- bind_rows(dfs)

write.csv(fulldf, "tblresults7-26-24_2.csv")