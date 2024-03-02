params.fastq = "fastq/*_{1,2}.fq.gz"
params.results = 'results'
params.aligner = 'bwamem2'
params.fastq_bwamem2 = false
params.fastq_bsbolt = false
params.min_mapq = 30
params.aligner_flags = '-SP5M'
params.assembly = 'hg38'
params.reference_genome_fasta = 'resources/index/bwamem2/hg38_noalts'
params.chromsizes = 'resources/hg38_noalts.sizes'
params.index_dir = 'resources/index/bwamem2'
params.index_prefix = 'hg38_noalts'
params.index_prefix_path = "${params.index_dir}/${params.index_prefix}"
params.all_index_files = []
params.index_ext_bwamem2 = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
params.select = """'(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")'"""

process INDEX_BWAMEM2 {
    tag "${sample_id}"
    publishDir params.index_dir, mode: 'copy'
    container "bskubi/bwa-mem2:latest"

    input:
    tuple val(index_prefix), path(reference_genome)

    output:
    tuple \
        val(index_prefix), \
        path("${index_prefix_path}.0123"), \
        path("${index_prefix_path}.amb"), \
        path("${index_prefix_path}.ann"), \
        path("${index_prefix_path}.bwt.2bit.64"), \
        path("${index_prefix_path}.pac")

    script:
    if (params.aligner == "bwamem2")
        """
        bwa-mem2 index ${reference_genome}
        """
}

process ALIGN {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'
    container "bskubi/bwa-mem2:latest"

    input:
    tuple val(sample_id), file(reads), path(index_dir), val(index_prefix)

    output:
    tuple val(sample_id), path("${sample_id}.bwamem2.bam")

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
    fastq_ch = Channel.fromFilePairs(params.fastq).map {sample_id, reads -> tuple(sample_id, reads, file(params.index_dir), params.index_prefix)}


    if (params.aligner == 'bwamem2') {
        index_exists = params.index_ext_bwamem2.collect {file(params.index_prefix_path + it).exists()}.every {it}
        reference_exists = file(params.reference_genome_fasta).exists()
        if (!index_exists) {
            print "bwa-mem2 index for ${params.assembly} not found, attempting to index based on ${params.reference_genome_fasta}"
            reference_genome_ch = Channel.fromPath(params.reference_genome_fasta, checkIfExists: true)
            index_input_ch = Channel.from(params.index_prefix).combine(reference_genome_ch)
            index_input_ch.view()
            INDEX_BWAMEM2(index_input_ch)
        }
    }
    
    fastq_ch.view()
    align_ch = ALIGN(fastq_ch)
    name_sort_ch = NAME_SORT(align_ch).map {sample_id, bam -> tuple(sample_id, bam, file(params.chromsizes))}
    PAIRS(name_sort_ch)
}
