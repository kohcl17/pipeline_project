#! /usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

df.list <- lapply(args, read.csv)
df.merge <- bind_rows(df.list)

write.csv(df.merge, "merged_sample_annotations.csv", row.names = FALSE)