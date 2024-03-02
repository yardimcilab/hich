params.fastq = "fastq/*_{1,2}.fq.gz"
params.results = 'results'
params.aligner = 'bwamem2'
params.fastq_bwamem2 = false
params.fastq_bsbolt = false
params.min_mapq = 30
params.aligner_flags = '-SP5M'
params.assembly = 'hg38'
params.index_prefix = 'hg38_noalts'
params.select = """'(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")'"""


process ALIGN {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'

    input:
    tuple val(sample_id), file(reads), path(index_dir), val(index_prefix)

    output:
    tuple val(sample_id), path("${sample_id}.bwamem2.bam")

    container "bskubi/bwa-mem2:latest"

    script:
    if (params.aligner == "bwamem2") 
        """
        bwa-mem2 mem ${params.aligner_flags} ${index_dir}/${index_prefix} ${reads[0]} ${reads[1]} | samtools view -b -o ${sample_id}.bwamem2.bam
        """
}

process NAME_SORT {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    container "clinicalgenomics/sambamba:0.8.0"
    
    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), file("${sample_id}.name_sort.bam")

    script:
    """
    sambamba sort -n -o ${sample_id}.name_sort.bam ${bam}
    """
}

process PAIRS {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    container "bskubi/pairtools:1.0.4"

    input:
    tuple \
        val(sample_id), \
        file(bam), \
        path(chromsizes)

    output:
    tuple \
        val(sample_id), \
        file("${sample_id}.selected.pairs.lz4"), \
        file("${sample_id}.not_selected.pairs.lz4"), \
        file("${sample_id}.parse.stats.txt"), \
        file("${sample_id}.select.stats.txt")

    script:
    """
    pairtools parse \
        ${bam} \
        -c ${chromsizes} \
        --assembly ${params.assembly} \
        --min-mapq ${params.min_mapq} \
        --output-stats ${sample_id}.parse.stats.txt \
    | \
    pairtools select \
        ${params.select} \
        --output-rest ${sample_id}.not_selected.pairs.lz4 \
    | \
    pairtools sort --output ${sample_id}.selected.pairs.lz4; \
    pairtools stats --output ${sample_id}.select.stats.txt ${sample_id}.selected.pairs.lz4
    """
}

workflow {
    fastq_ch = Channel.fromFilePairs(params.fastq).map {sample_id, reads -> tuple(sample_id, reads, params.index_dir, params.index_prefix)}
    align_ch = ALIGN_BWAMEM2(fastq_ch)
    name_sort_ch = NAME_SORT(align_ch).map {sample_id, bam -> tuple(sample_id, bam, params.chromsizes)}
    PAIRS(name_sort_ch)
}
