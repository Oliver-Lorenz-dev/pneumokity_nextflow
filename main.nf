#!/usr/bin/env nextflow

process PNEUMOKITY {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.csv"
    tag "$sample_id"
    container '/software/pathogen/images/pneumokity-1.0.simg'
    input:
    tuple val(sample_id), path(read_1), path(read_2)

    output:
    path("*.csv"), emit: results_ch
    tuple val(sample_id), path("workdir"), emit: workdir_ch
    script:
    """
    pwd > workdir
    pneumokity.py mix -f ${read_1} ${read_2} -m /opt/conda/bin/mash -o .
    mv pneumo_capsular_typing/*_result_data.csv .
    """
}

process GET_RUN_TIME {
    publishDir "${params.results_dir}", mode: 'copy', overwrite: true, pattern: "*.txt"
    tag "$sample_id"
    input:
    tuple val(sample_id), path(workdir)

    output:
    path("*.txt"), emit: runtime_ch
    script:
    """
    pneumokity_workdir=\$(cat ${workdir})
    grep "Run time" \${pneumokity_workdir}/.command.log | awk '{ print \$(NF-1) }' > ${sample_id}_run_time.txt
    """
}

workflow {
    manifest_ch = Channel.fromPath(params.manifest)

    fastq_path_ch = manifest_ch.splitCsv(header: true, sep: ',')
            .map{ row -> tuple(row.sample_id, file(row.read_1), file(row.read_2)) }

    PNEUMOKITY(fastq_path_ch)

    GET_RUN_TIME(PNEUMOKITY.out.workdir_ch)
}