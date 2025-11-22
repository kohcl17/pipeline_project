#! /usr/bin/env Rscript

library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- file.path(".", "filtered_contigs.csv")

filter.by <- str_replace_all(args[2], ",", "|")

df <- read.csv(infile)
x2 <- filter(df, grepl(x = c_gene, pattern = filter.by))

write.csv(x2, outfile, row.names = FALSE)