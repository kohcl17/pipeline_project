// run pyir
process run_pyIR {
    // this only publishes the directory AFTER the process
    publishDir "${params.output}", mode: 'copy', overwrite: true

    // take in the input file paths and the filters to use
    input:
    tuple val(sampleID), path(annotation), path(input_fasta)

    output:
    tuple val(sampleID), path(annotation), path("${sampleID}_pyir_output.tsv")

    script:
    """
    pyir ${input_fasta} \
        --enable_filter \
        -o ${sampleID}_pyir_output \
        --outfmt tsv \
        -z gunzip
    """
}

// merged pyIR output and adds the V3J clonotypes
process merge_and_annotate {
    // this only publishes the directory AFTER the process
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    // take in the input file paths and the filters to use
    input:
    tuple val (sampleID), path(annotation), path(pyirOut)

    // push the csv file. If just using * it will also recursively push the directory
    output:
    path "*.csv"

    // the R script writes INTO working directory (ie. when Nextflow starts it says 'work')
    script:
    """
    merge_and_annotate.R ${sampleID} ${annotation} ${pyirOut}
    """
}

process merge_samples {
    // this only publishes the directory AFTER the process
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    // take in all the csvs from merge_and_annotate
    input:
    path(csvs)

    // push the csv file
    output:
    path "merged_sample_annotations.csv"


    script:
    """
    merge_samples.R ${csvs}
    """
}

// filter c_gene column (heavy chain)
process filter_IGH_contigs {
    // this only publishes the directory AFTER the process
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    // take in the input file paths and the filters to use
    input:
    path input_file
    val filterBy

    // push the csv file. If just using * it will also recursively push the directory
    output:
    path "*.csv"

    script:
    """
    filter_IGH_contigs.R ${input_file} ${filterBy}
    """
}

// TODO: think about whether to ask for unique clonotypes or common?
// process find_unique_clonotypes {
//     // this only publishes the directory AFTER the process
//     publishDir "${params.output}", mode: 'copy', overwrite: true
    
//     // take in the input file paths and the filters to use
//     input:
//     path input_file

//     // push the csv file. If just using * it will also recursively push the directory
//     output:
//     path "*.csv"

//     script:
//     """
//     find_unique_clonotypes.R ${input_file}
//     """
// }

// Do multiple sequence alignment
process msa {
    // this only publishes the directory AFTER the process
    publishDir "${params.output}", mode: 'copy', overwrite: true
    
    // take in the input file paths and the filters to use
    input:
    path input_file
    val alignWhich
    val doPCOA

    // push the rds file. If just using * it will also recursively push the directory
    output:
    path "*.rds"
    path "*.pdf", optional: true

    script:
    """
    msa.R ${input_file} ${alignWhich} ${doPCOA}
    """
}

// Workflow block
workflow {
    ch_map_file = channel.\
    fromPath(params.inputMap).\
    splitCsv(sep: ",").\
    map{ row -> [row[0], file(row[1]), file(row[2]) ] }

    // nFiles = ch_map_file.count()
    
    pyirOut = run_pyIR(ch_map_file)
    
    merged_file = merge_and_annotate(pyirOut)
    ch_merged_files = merged_file.collect()
    merge_sample_files = merge_samples(ch_merged_files)

    // Manage the filteredContigs argument
    if (params.filterIGHContigs instanceof String) {
        filterChoice = params.filterIGHContigs.toUpperCase()
        filtered_file = filter_IGH_contigs(merge_sample_files, filterChoice)

    } else {
        filtered_file = merge_sample_files
    }
    
    // if (nFiles > 1) {
    //     unique_clonotype_files = find_unique_clonotypes(filtered_file)
    // } else {
    //     unique_clonotype_files = filtered_file
    // }

    msaOut = msa(filtered_file, params.alignWhich, params.doPCoA)
}