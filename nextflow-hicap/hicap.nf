#!/usr/bin/env nextflow
/*
 * Define the default parameters 
 */
params.data = "${baseDir}/data"
params.sample_ids = "${baseDir}/data/sample_ids"
params.outdir = "${baseDir}/results"
params.primer_info = "${baseDir}/data/primer_info"
params.script_directory = "${baseDir}/bin"

log.info """ 
    H I C A P - T Y P I N G - P I P E L I N E
    ========================================
    Data directory      : ${params.data}
    Sample IDs          : ${params.sample_ids}
    Output directory    : ${params.outdir}
    Primer info         : ${params.primer_info}
    Script directory    : ${params.script_directory}
    """
    .stripIndent()

/*
 * The `locate_fasta` process generates the path for each sample ID using `mdu contigs -s sample_ID`.
 */
process locate_fasta {
    input:
    path sampleIDs

    output:
    path "fasta_path.txt"

    script:
    """
    while read sample; do
        # Get the path to the FASTA file
        mdu contigs --sample_id "\$sample" >> fasta_path.txt
    done < ${sampleIDs}
    """
}

/*
 * The `copy_fasta` process copies the FASTA file to the `fasta_files` directory within the data folder with the sample ID name.
 */

process copy_fasta_run_hicap {
    publishDir path: {params.outdir}, mode: 'copy'
    // beforeScript 'source activate /home/himals/.conda/envs/hicap_dev'
    conda 'env.yml'
    
    input:
    path fasta_path

    output:
    path "output/"
    file "hicap_summary.csv"

    shell:
    '''
    mkdir -p fasta_files/

    while read path; do
       # Assign the sample ID and the path to the FASTA file to variables
       sample_id=$(echo $path | awk '{print $1}')
       fasta_path=$(echo $path | awk '{print $2}')
       
       # Copy the FASTA file to the fasta_files directory within the data folder with the sample ID name
       cp $fasta_path fasta_files/$sample_id.fasta     
    done < !{fasta_path}

    mkdir -p output

    # Execute the command for each FASTA file
    parallel hicap --query_fp {} --output_dir output/ ::: fasta_files/*fasta

    # Generate .tsv files for samples that didn't have any hits
    # This is necessary because the hicap command doesn't generate a .tsv file if there are no hits
    # This will cause the hicap_summary.tab file to be missing some samples
    
    while read path; do
        sample_id=$(echo $path | awk '{print $1}')

        file_path="output/${sample_id}.tsv"
        if [[ ! -f "$file_path" ]]; then
            # The file doesn't exist, generate sample_id.tsv
            echo -e "#isolate\tpredicted_serotype\tattributes\tgenes_identified\tlocus_location\tregion_I_genes\tregion_II_genes\tregion_III_genes\tIS1016_hits" > "$file_path"
            echo -e "${sample_id}\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found\tno hits to any cap locus gene found" >> "$file_path"
        fi
    done < !{fasta_path}

    # Combine the output
    summary_fps=(output/*tsv);
    { head -n1 "${summary_fps[0]}"; sed '/^#/d' "${summary_fps[@]}"; } > hicap_summary.tab
    
    # Convert tsv to csv
    sed '1s/^#//' hicap_summary.tab | csvtk tab2csv > hicap_summary.csv
    '''
}

/*
 * the `seqkit_pcr` process uses the `seqkit pcr` command to extract the region of the FASTA file that is between the primer sequences.
 */

process seqkit_pcr {
    publishDir path: {params.outdir}, mode: 'copy'
    conda 'env.yml'
    
    input:
    path sampleIDs
    path primer_info

    output:
    path "detection_results.csv"

    shell:
    '''
    seqkit-pcr.sh !{primer_info} !{sampleIDs} 
    '''
}

process spread_data {
    publishDir path: {params.outdir}, mode: 'copy'
    conda 'env.yml'

    input:
    path script_directory
    path detection_results

    output:
    path "detection_results_wide.csv"

    shell:
    '''

    python !{script_directory}/spread.py !{detection_results} detection_results_wide.csv
    '''
}

workflow {
    locate_fasta(params.sample_ids)
    copy_fasta_run_hicap(locate_fasta.out)
    seqkit_pcr(params.sample_ids, params.primer_info)
    spread_data(params.script_directory, seqkit_pcr.out)
}




