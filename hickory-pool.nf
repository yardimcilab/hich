include { CAT } from './hickory-common.nf'
params.results = 'results'

process DEDUP {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    container "bskubi/pairtools:1.0.4"

    input:
    tuple \
        val(block_id), \
        val(sample_id), \
        path(pairs),
        val(timing)
    
    output:
    tuple \
        val(block_id), \
        val(sample_id), \
        path("${sample_id}.dedup_${timing}.pairs.lz4"), \
        path("${sample_id}.dups_${timing}.pairs.lz4"), \
        path("${sample_id}.dedup_${timing}.stats.txt")

    script:
    """
    pairtools dedup \
        --send-header-to both \
        --mark-dups \
        --output ${sample_id}.dedup_${timing}.pairs.lz4 \
        --output-dups ${sample_id}.dups_${timing}.pairs.lz4 \
        --output-stats ${sample_id}.dedup_${timing}.stats.txt \
        ${pairs}
    """
}

process MERGE {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    container "bskubi/pairtools:1.0.4"

    input:
    tuple \
        val(block_id), \
        val(sample_ids), \
        path(files)
    
    output:
    tuple \
        val(block_id), \
        path("${block_id}.merge.pairs.lz4")

    script:
    """
    pairtools merge --output "${block_id}.merge.pairs.lz4" ${files}
    """
}

workflow final_dedup {
    take: merge_ch

    main:
    DEDUP(merge_ch.map{tuple(it[0], it[0], it[1], "after")})
}

workflow {
    /*
    The workflow needs to know
        What the input .pairs paths are
        How to associate .pairs paths with sample ids.
        What the block IDs are
        How to associate sample IDs with block IDs

        We can plausibly associate .pairs paths with sample ids and sample IDs block IDs in any of the following ways:
            - A database with a specific schema from which we obtain a join table
            - A join table (from CSV or a database) with FILENAME, SAMPLE_ID, and BLOCK_ID entries
            - A regex extracting SAMPLE_ID from FILENAME and BLOCK_ID from SAMPLE_ID
            - JSON that stores these associations (i.e. key = BLOCK_ID, value = list of (FILENAME, SAMPLE_ID) pairs)
        
        To extract from CSV, we need to read each row and store it in a channel
    */
    csv_ch = CAT(Channel.fromPath('demo.csv')).splitCsv(skip: 1).map {tuple(it[0], it[1], file(it[2]), "before")}
    block_id_ch = csv_ch.map{it[0]}.distinct()

    dedup_ch = DEDUP(csv_ch).map{tuple(it[0], it[1], it[2])}.groupTuple()
    merge_ch = MERGE(dedup_ch)
    final_dedup(merge_ch)
}