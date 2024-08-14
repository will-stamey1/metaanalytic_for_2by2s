library(dplyr)

matching_files <- list.files(pattern = "tbl1and2sims_")

print(matching_files)

dfs <- lapply(matching_files, read.csv)

#for(f in 1:length(matching_files)){
#  print(paste(matching_files[f], dim(dfs[[f]])))
#} 

fulldf <- bind_rows(dfs)

write.csv(fulldf, "tblresults8-14-24_unifsdprior.csv")