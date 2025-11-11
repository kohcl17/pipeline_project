# Dataset
The dataset was from [CellRanger vdj example dataset](https://www.10xgenomics.com/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0). The dataset contains 1105 CD19+ B-cells from a healthy male donor aged 27. The filtered contig annotations csv file contains the annotated VDJ contigs that passed cellranger's filtering algorithm. The filtered fasta files should contain the same number of annotations as the csv file but in fasta format (ie. just sequences and the contig sequence id).

# Environment

## Python
python==3.13.5\
pyir==1.5.0

## R
R==4.5.1\
tidytree==0.4.6\
Biostrings==2.76.0\
tibble==3.3.0\
msa==1.40.0\
ggplot2==3.5.2\
stringr==1.5.2\
ggtree=3.16.3\
RColorBrewer==1.1-3\
dplyr==1.1.4

# Pipeline Description
This pipeline is meant to take as input the filtered fasta and contig annotation files from cellranger's vdj pipeline, and will output the multiple sequence alignment R data along with optional principal coordinate analysis (PCoA) plots of the sequence alignments. The pipeline will first annotate the contig annotation dataset with V3J clonotype information, denoted by the combination of a V gene, J gene, and CDR3 amino acid sequence, as described in this article by [Soto et al (2019)](https://www.nature.com/articles/s41586-019-0934-8).

# Usage

To the `nextflow.config` file, use the following configuration options:
```
params {
    inputMap = './input.csv'
    output = "./results"
    filterIGHContigs = false
    alignWhich = "DNA"
    doPCoA = true
}
```
Then use the pipeline by entering the following into the terminal:

```
nextflow run main.nf
```

Another way to specify these options is to directly call it in the command line:

```
nextflow run main.nf \
    --inputMap ./input.csv \
    --output ./results \
    --filterIGHContigs false \
    --alignWhich DNA \
    --doPCoA true
```

# TODO: add the conda env or docker/sif image info