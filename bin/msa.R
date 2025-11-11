#! /usr/bin/env Rscript

# basic utilities
library(dplyr)
library(RColorBrewer)
library(ggplot2)

# Multiple sequence alignment
library(Biostrings)
library(ggtree)
library(msa)
library(tidytree)


############################## FUNCTIONS ########################################

# Working assumptions
## 1. Everything in the sequence.df will be aligned so do the filtering first
## 2. Standard identity matrix will be used for distance alignment
## 3. ape's neighbour joining tree will be used for tree creation

# Function params
# sequence.df: the merged IGHG dataset filtered for the samples to be aligned
# col.alignment: string of column in sequence.df to use for alignment, eg. 'sequence_alignment'
# seq.type: 'DNA' or 'AA'
# col.seq.name: string of column to use for names of sequences (eg. sequence_id)
align.sequences <- function(sequence.df, 
                            col.alignment, 
                            seq.type = "DNA", 
                            col.seq.name = "sequence_id"
                            ) {
  rename.ls <- c("label" = col.seq.name)
  sequence.df <- rename(sequence.df, any_of(rename.ls))
  
  seq.vct <- sequence.df[[col.alignment]]
  names(seq.vct) <- sequence.df[['label']]
  
  if (seq.type == "DNA") {
    seq.strSet <- DNAStringSet(seq.vct)
  } else if (seq.type == "AA") {
    seq.strSet <- AAStringSet(seq.vct)
  }
  
  seq.msa <- msa(seq.strSet, method = "ClustalOmega")
  seq.msa <- msaConvert(seq.msa, type="seqinr::alignment")
  
  dist.mtx <- seqinr::dist.alignment(seq.msa, matrix = "identity")
  
  seq.tree.ape <- ape::nj(dist.mtx)
  
  seq.tree <- full_join(as.treedata(seq.tree.ape),
                        sequence.df,
                        by = "label"
                        )
  return(list(tree = seq.tree,
              tree.ape = seq.tree.ape,
              msa = seq.msa, 
              dist.mtx = dist.mtx, 
              strSet = seq.strSet,
              seqDf = sequence.df))
}

PCOA_plot <- function(seq.file) {
  seqDf <- seq.file$seqDf
	tree <- seq.file$tree
  dist.data <- as.matrix(seq.file$dist.mtx)
  tree.coords <- ggtree(tree)$data %>% 
                  filter(isTip == TRUE) %>%
                  select(label,x,y)

  pcoa <- cmdscale(dist.data, eig = TRUE, k = 2)
  xy.pcoa <- data.frame(pcoa$points, check.names = FALSE)
  colnames(xy.pcoa) <- c('PCoA1', 'PCoA2')
  
  combined.data <- merge(xy.pcoa,
                          seqDf, 
                          by.x = 0, 
                          by.y = 'label',
                          all.x = TRUE,
                          all.y = FALSE)
                  
  combined.data <- merge(combined.data, tree.coords, by.x = "Row.names", by.y = "label")
  
  pcoa.plot.clr.dist <- ggplot(combined.data,
                               aes(x = PCoA1, 
                                   y = PCoA2, 
                                   colour = x)
                      ) +
    geom_point() +
    labs(title = nm) +
    scale_colour_gradientn(colours = brewer.pal(6, "YlGnBu"),
                          name = "distance from 0"
                          ) +
    theme_classic()

    return(pcoa.plot.clr.dist)
}

#################################################################################
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
alignWhich <- args[2]
doPCOA <- args[3]

col.alignment <- ifelse(alignWhich == "DNA",
    "sequence_alignment",
    "sequence_alignment_aa"
)

outfile <- file.path(".", "alignment_data.rds")

df.list <- read.csv(infile)
df.list <- split(df.list, ~nextflow.sample.id)
df.list <- lapply(df.list, ungroup)

aligned.list <- list()

for (nm in names(df.list)) {
	df <- df.list[[nm]]
	seq.data <- align.sequences(df,
		col.alignment = col.alignment,
		col.seq.name = "sequence_id",
		seq.type = alignWhich
	)
	if (doPCOA) {
    plot.name <- paste("./", nm, "_PCoA.pdf", sep = "")
		pcoa.plot <- PCOA_plot(seq.data)
    ggsave(plot.name, pcoa.plot, width = 5, height = 4)
	}
 aligned.list[[nm]] <- seq.data
}

saveRDS(aligned.list, outfile)