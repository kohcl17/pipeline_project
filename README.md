# Dataset
The dataset was from the [CellRanger vdj example dataset](https://www.10xgenomics.com/datasets/human-b-cells-from-a-healthy-donor-1-k-cells-2-standard-6-0-0). The dataset contains 1105 CD19+ B-cells from a healthy male donor aged 27. The filtered contig annotations csv file contains the annotated VDJ contigs that passed cellranger's filtering algorithm. The filtered fasta files should contain the same number of annotations as the csv file but in fasta format (ie. just sequences and the contig sequence id).

# Hypothesis

I hypothesise that majority of the aligned sequences that are filtered by IGHG constant region genes will predominantly be composed of the IGHG1 subgroup, as this subgroup predominantly targets toxins and viruses. Seeing as the dataset is from a healthy individual, I expect to see a modest B-cell receptor (BCR) repertoire as less clonal expansion would have occurred, thus resulting in smaller BCR contig mutation distances from the phylogenetic tree and lower diversity in terms of VDJ profiles.

# Aim 
I aim to identify the most mutated contig sequences with the highest number of unique molecular identifiers (UMIs) as this indicates higher diversity and possibly clonal expansion. 

# Objective
The objective of this pipeline is to identify BCR sequences that contribute to greater clonal diversity in the hopes that running it with a particular disease perturbation, for example H5N1, will yield a set of sequences that can be used for monoclonal antibody production.

# Pipeline Description

The following steps will be taken in this pipeline:

1. Annotate the contig annotation dataset with V3J clonotype information, denoted by the combination of a V gene, J gene, and CDR3 amino acid sequence, as described in this article by [Soto et al (2019)](https://www.nature.com/articles/s41586-019-0934-8). 

2. [PyIR](https://doi.org/10.1186/s12859-020-03649-5), a Python wrapper for IgBLAST, will be used to filter out lower quality annotations such as contigs that contain stop codons, CDR3 regions that are non-continuous, and junctional regions that are out of frame [Gilchuk et al (2020)](https://www.nature.com/articles/s41551-020-0594-x#Sec10). 

3. (Optional) Filter the contig annotations based on user specifications.

4. Merge and annotate the PyIR output with the cellranger annotations.

5. Perform multiple sequence alignment (MSA) with the ClustalOmega algorithm.

6. Graphs of the phylogenetic tree alignments will be generated along with an optional principal coordinate analysis (PCoA) plot for each sample.

7. The filtered contigs will be ranked by mutation distance as determined by MSA as well as number of clones based on UMIs.

The graph below shows where the configuration parameters will be used and the sequence of processes in this pipeline.



## Input

input.csv file with the following specifications:
- sampleName, absolute path to filtered fasta, absolute path to filtered contig annotation files

> [!NOTE]
> Both fasta and contig annotation files should be provided by CellRanger output from its vdj pipeline.

## Output 

**pyir/**
- sample_pyir_output.tsv: file containing pyIR output

**merge_and_annotate/**
- each sample's pyIR output merged with cellranger contig annotation files
- this is prior to filtering the samples by c_gene column, an optional parameter to be set

**filtered/**
- sample-merged filtered contigs based on the optional parameter --filterIGHContigs

**msa/**
- alignment_data.rds: an R data structure containing a list of all the samples' multiple sequence alignment (MSA) data and supporting files
- sample_MSA_tree.pdf: phylogenetic tree of MSA
- sample_PCoA.pdf: Principal Coordinate Analysis (PCoA) plots of the MSA tree as an alternative way to visualise the data

**clonotype_rank/**
- sample_rank_paired.csv: a file containing an additional column 'rank' which ranks the contigs

pipeline_dag.html: the directed acrylic graph (DAG) showing how the pipeline has used the parameters and input files, as well as how the processes link to each other.

# Environment

## Nextflow/Java
openjdk==17.0.10\
nextflow==25.04.8.5956

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

# Usage

For compatibility with most computers, once this repo has been downloaded or cloned, build a docker image with the following command in the terminal:

> [!NOTE]
> Please ensure that docker is running in the background before running the pipeline.

```
cd path/to/repo
docker build -t pipeline-501a:latest . --platform linux/amd64
```

If you run into errors regarding permissions, run the following lines in the terminal:

```
chmod +x bin/*
```

To the `nextflow.config` file, use the following configuration options:
```
params {
    inputMap = './input.csv'
    output = "./results"
    filterIGHContigs = 'IGHG'
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
    --filterIGHContigs IGHG \
    --alignWhich DNA \
    --doPCoA true \
    --colourMSATipsBy c_gene \
    -profile docker
```

# Parameters
## Input/Output
**inputMap**: str | path to input csv file

**output**: str | path to create an output directory

## Contig Filtering & Alignment
**filterIGContigs**: false or str regex pattern of IG contigs in the c_gene column of filtered_contig_annotations.csv | filters the c_gene column

**alignWhich**: str | DNA or AA

## Graphs
**doPCoA**: bool | add PCoA graphs which allows for condensed visualisation of the MSA phylogenetic tree

**colourMSATipsBy**: str | colour the MSA phylogenetic tree by a certain column. During the annotation part of the pipeline, a IGHV subgroup column will be created and can be used for colouring by specifying IGHV.subgrp

## Nextflow parameters
**profile**: if no profile has been specified, nextflow will use your local environment which may not have the necessary packages installed. For best performance, create a docker image from the instructions above and specify profile docker.

# Expected Output
The expected output can be found in the expected_results folder of this repository.

# Future Work

- Basic visualisation of repertoire diversity, with some emphasis on V3J clonotypes
- Contigs from multiple samples or even group-level comparisons