#! /usr/bin/env Rscript

# basic utilities
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)

# Multiple sequence alignment
library(Biostrings)
library(ggtree)
library(msa)
library(tidytree)

############ FUNCTIONS #############
get.xdist <- function(treedat) {
  dist.dat <- ggtree(treedat)$data %>%
    filter(isTip == TRUE) %>%
    select(label, x, y)
  return(dist.dat)
}

rank.and.pair <- function(dist.dat, filtered.annot, merged.sample.annotations) {
  cdat <- merge(dist.dat,
                filtered.annot,
                by.x = 'label', 
                by.y = 'sequence_id',
                all.x = TRUE,
                all.y = FALSE) %>%
    arrange(desc(x))

  cdat.rank <- cdat %>%
    arrange(desc(x), desc(Clones), .by_group = TRUE) %>%
    mutate(rank = 1:nrow(.)) %>%
    select(barcode, rank)

  ranked.pairs <- merged.sample.annotations %>%
    filter(barcode %in% cdat.rank$barcode) %>% # only use the pairs that were in the filtered annotations
    left_join(y = cdat.rank, by = 'barcode') %>%
    arrange(rank)
  return(ranked.pairs)
}


####################################
args <- commandArgs(trailingOnly = TRUE)
infile.rds <- args[1]
infile.sample.annot <- args[2]
infile.filtered.annot <- args[3]

# READ DATA
seq.data <- readRDS(infile.rds)

sample.annot <- read.csv(infile.sample.annot)
sample.annot.list <- split(sample.annot, ~nextflow.sample.id)
sample.annot.list <- lapply(sample.annot.list, ungroup)

filtered.annot <- read.csv(infile.filtered.annot)
filtered.annot.list <- split(filtered.annot, ~nextflow.sample.id)
filtered.annot.list <- lapply(filtered.annot.list, ungroup)

# GET DISTANCES
dist.dats <- lapply(seq.data, function(x) {
  tr_obj <- x$tree
  return(get.xdist(tr_obj))
  })

# RANK AND PAIR ANNOTATIONS
for (nm in names(dist.dats)) {
  dst <- dist.dats[[nm]]
  annot <- sample.annot.list[[nm]]
  f.annot <- filtered.annot.list[[nm]]
  out.name <- file.path(".", paste(nm, "rank_paired.csv", sep = "_"))
  rnk <- rank.and.pair(dst, f.annot, annot)
  write.csv(rnk, out.name, row.names = FALSE)
}


# debugging purposes
# write.xlsx(ranked.seqs, "./out_ranked_seqs.xlsx")