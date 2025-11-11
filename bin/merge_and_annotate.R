#! /usr/bin/env Rscript

library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
sample.id <- args[1]
annot.file <- args[2]
pyir.file <- args[3]

outfile <- file.path(".", paste(sample.id, "v3j_annotations.csv", sep = "_"))

annot.df <- read.csv(annot.file)
annot.df <- mutate(annot.df, v3j = paste(v_gene, j_gene, cdr3, sep = "_")) # Note here this is the cellranger format

pyir.df <- read.table(pyir.file, sep = "\t", header = TRUE)

annot.df <- annot.df |>
    dplyr::select(barcode,
                 contig_id,
                 v3j,
                 c_gene, 
                 reads,
                 umis, 
                 raw_clonotype_id,
                 raw_consensus_id,
                 exact_subclonotype_id
                 )

df.merge <- merge(pyir.df, 
                  annot.df,
                  by.x = 'sequence_id',
                  by.y = 'contig_id',
                  all.x = TRUE,
                  all.y = FALSE
                )

df.merge <- df.merge |>
    rename(Clones = umis,
           CDR3.nt = cdr3,
           CDR3.aa = cdr3_aa,
           V.name = v_family,
           D.name = d_family,
           J.name = j_family,
           V.end = v_alignment_end,
           D.start = d_alignment_start,
           D.end = d_alignment_end,
           J.start = j_alignment_start,
           J.end = j_alignment_end
        )
df.merge <- df.merge |>
    mutate(nextflow.sample.id = sample.id, 
          IGHV.subgrp = str_replace_all(V.name, "-.*", ""))

write.csv(df.merge, outfile, row.names = FALSE)