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
dplyr==1.1.4 \
patchwork==1.3.2

# Pipeline Description
This pipeline takes as input the filtered fasta and contig annotation files from cellranger's vdj pipeline, and will output the multiple sequence alignment (MSA) R data along with the phylogenetic tree and an optional principal coordinate analysis (PCoA) plots of the sequence alignments. The pipeline will first annotate the contig annotation dataset with V3J clonotype information, denoted by the combination of a V gene, J gene, and CDR3 amino acid sequence, as described in this article by [Soto et al (2019)](https://www.nature.com/articles/s41586-019-0934-8). Next, PyIR will be used to filter out lower quality annotations such as contigs that contain stop codons, CDR3 regions that are non-continuous, and junctional regions that are out of frame [Gilchuk et al (2020)](https://www.nature.com/articles/s41551-020-0594-x#Sec10). Thereafter, the pipeline will optionally filter the contig annotations based on user specifications, merge and annotate the PyIR output with the cellranger annotations, and perform MSA with the ClustalOmega algorithm.

# Usage

For compatibility with most computers, once this repo has been downloaded or cloned, build a docker image with the following command in the terminal:

> [!NOTE]
> Please ensure that docker is running in the background before running the pipeline.

```
cd path/to/repo
docker build -t pipeline-501a:latest . --platform linux/amd64
```

To the `nextflow.config` file, use the following configuration options:
```
params {
    inputMap = './input.csv'
    output = "./results"
    filterIGHContigs = false
    alignWhich = "DNA"
    doPCoA = true
    colourMSATipsBy = 'c_gene'
}
```
Then use the pipeline by entering the following into the terminal:

```
nextflow run main.nf -profile docker
```

Another way to specify these options is to directly call them in the command line:

```
nextflow run main.nf \
    --inputMap ./input.csv \
    --output ./results \
    --filterIGHContigs false \
    --alignWhich DNA \
    --doPCoA true \
    --colourMSATipsBy c_gene \
    -profile docker
```

# Parameters
## Input/Output
**inputMap**: str | path to csv file containing sample_name,path/to/cellranger/filtered_contig_annotations.csv,path/to/cellranger/filtered_contig.fasta

**output**: str | path to create an output directory

## Contig Filtering & Alignment
**filterIGContigs**: false or str of IG contigs in the c_gene column of filtered_contig_annotations.csv | filters the c_gene column

**alignWhich**: str | DNA or AA

## Graphs
**doPCoA**: bool | add PCoA graphs which allows for condensed visualisation of the MSA phylogenetic tree

**colourMSATipsBy**: str | colour the MSA phylogenetic tree by a certain column. During the annotation part of the pipeline, a IGHV subgroup column will be created and can be used for colouring by specifying IGHV.subgrp

## Nextflow parameters
**profile**: if no profile has been specified, nextflow will use your local environment which may not have the necessary packages installed. For best performance, create a docker image from the instructions above and specify profile docker.
