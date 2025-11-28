#! /usr/bin/env Rscript

library(dplyr)
library(stringr)
library(tibble)

input.files <- commandArgs(trailingOnly = TRUE)[1]

df <- read.csv(input.files)
df.list <- df |>
    group_split(nextflow.sample.id)

all.comb.mat <- make_comb_mat(lapply(df.list, function(x) x$v3j))
all.upset <- UpSet(all.comb.mat,
                   column_title = "Unique Clonotypes",
                   top_annotation = upset_top_annotation(all.comb.mat, 
                                                         add_numbers = TRUE,
                                                         numbers_rot = 0),
                   right_annotation = upset_right_annotation(all.comb.mat,
                                                              add_numbers = TRUE)
                   )