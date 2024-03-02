params.results = 'results'

process MULTIQC {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    conda 'multiqc.yaml'

    input:
    tuple \
        val(partition_id), \
        path(stats)
    
    output:
    path("${partition_id}_multiqc_report.html")

    script:
    """
    multiqc -m pairtools -f -o ${partition_id}_multiqc_report.html ${stats}
    """
}

workflow {
    stats_ch = Channel.fromPath("${params.results}/*stats.txt")
    multiqc_ch = Channel.from("all").combine(stats_ch).groupTuple()

    MULTIQC(multiqc_ch)
}