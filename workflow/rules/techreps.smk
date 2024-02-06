rule bwamem2_bam:
    """
    Align paired-end reads using BWA MEM2 with parameters from 4D Nucleome Hi-C Processing Pipeline guidelines
    Then convert to .bam using samtools.

    This aligns each end singly, calling the 5'-most alignments as primary and annotates
    secondary/supplementary clipped reads as secondary.

    With this approach, it is not necessary to do junction trimming, and in fact this is not advised.
    Junction trimming will artificially make complex walks appear as U/U reads. Complex walks that
    are 'rescued' in this way typically appear to be the result of assay noise. They are excluded by
    the pipeline.

    See https://pairtools.readthedocs.io/en/latest/parsing.html#rescuing-complex-walks

    As alignment speed scales almost linearly with number of processes, we align one file at a time
    to save on memory.
    https://github.com/bwa-mem2/bwa-mem2
    https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline
    https://github.com/samtools/samtools
    """
    input:
        r1 = config['data']['fastq'] + "{cond}_{br}_{tr}_1.fq.gz",
        r2 = config['data']['fastq'] + "{cond}_{br}_{tr}_2.fq.gz",
        index = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['index']
    output: "results/{proc}/{cond}/{br}/{tr}/{tr}.bam"
    params:
        min_mapq = lambda wc: config['processes'][wc.proc]['parameters']['min_mapq']['align']
    container:
        "docker://bskubi/bwa-mem2:latest"
    benchmark: "benchmark/bwamem2_bam/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/bwamem2_bam/{proc}/{cond}/{br}/{tr}.log"
    threads: workflow.cores
    shell:
        """
        bwa-mem2 mem -t {threads} \
        -SP5M \
        -T {params.min_mapq} \
        -P {input.index} \
        {input.r1:q} {input.r2:q} \
        2>> {log} \
        | samtools view -b -o {output:q} 2>> {log}
        """

rule name_sort:
    """
    Sort .bam files by name with sambamba
    https://github.com/biod/sambamba
    """
    input: "results/{proc}/{cond}/{br}/{tr}/{tr}.bam"
    output: "results/{proc}/{cond}/{br}/{tr}/{tr}_name_sort.bam"
    resources:
        mem = "2GiB"
    container: "docker://clinicalgenomics/sambamba:0.8.0"
    benchmark: "benchmark/name_sort/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/name_sort/{proc}/{cond}/{br}/{tr}.log"
    shell: """sambamba sort -n -m 2GiB -t {threads} -o {output:q} {input:q} &>> {log}"""

rule parse:
    """
    Convert .bam to .pairs with pairtools. See https://github.com/open2c/pairtools

    Also outputs statistics file for inclusion in multiqc report.
    """
    input:
        bam = "results/{proc}/{cond}/{br}/{tr}/{tr}_name_sort.bam",
        chromsizes = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['chromsizes']
    output:
        pairtools_parse = "results/{proc}/{cond}/{br}/{tr}/{tr}_parse.pairs",
        pairtools_parse_stats = "results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_parse_stats.txt"
    params:
        assembly = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly'],
        min_mapq = lambda wc: config['processes'][wc.proc]['parameters']['min_mapq']['parse']
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/parse/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/parse/{proc}/{cond}/{br}/{tr}.log"
    threads: workflow.cores
    shell:
        """
        pairtools parse {input.bam:q} \
            --min-mapq {params.min_mapq} \
            -o {output.pairtools_parse:q} \
            --output-stats {output.pairtools_parse_stats:q} \
            -c {input.chromsizes:q} \
            --assembly {params.assembly} \
            --nproc-in {threads} \
            --nproc-out {threads} \
            &>> {log}
        """

rule select:
    """
    Select only valid pairs (with 2x alignments called as unique ('U') or rescued by pairtools ('R'))
    for single junction ligations. Exclude complex walks.
    """
    input: "results/{proc}/{cond}/{br}/{tr}/{tr}_parse.pairs"
    output:
        select = "results/{proc}/{cond}/{br}/{tr}/{tr}_select.pairs",
        rest = "results/{proc}/{cond}/{br}/{tr}/{tr}_rest.pairs"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/select/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/select/{proc}/{cond}/{br}/{tr}.log"
    shell:
        """pairtools select '(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")' \
                    --output-rest {output.rest:q} \
                    -o {output.select:q} \
                    {input:q} \
                    &>> {log}"""

rule pos_sort_tr:
    """
    Sort .pairs by position
    """
    input: "results/{proc}/{cond}/{br}/{tr}/{tr}_select.pairs"
    output: "results/{proc}/{cond}/{br}/{tr}/{tr}_pos_sort.pairs"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/pos_sort/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/pos_sort/{proc}/{cond}/{br}/{tr}.log"
    shell: """pairtools sort -o {output:q} {input:q} &>> {log}"""

rule dedup_tr:
    """
    Deduplicate technical replicates *prior* to merging (biological replicates
    wil be deduplicated *after* merging). Also outputs statistics file for
    inclusion in multiqc report.
    """
    input: "results/{proc}/{cond}/{br}/{tr}/{tr}_pos_sort.pairs"
    output:
        dedup = "results/{proc}/{cond}/{br}/{tr}/{tr}_dedup.pairs",
        dups = "results/{proc}/{cond}/{br}/{tr}/{tr}_duplicates.pairs",
        stats = "results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_dedup_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/dedup/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/dedup/{proc}/{cond}/{br}/{tr}.log"
    shell:
        """
        pairtools dedup \
            --send-header-to both \
            --mark-dups \
            --output-stats {output.stats:q} \
            --mark-dups \
            --output-dups {output.dups:q} \
            -o {output.dedup:q} \
            {input:q} \
            &>> {log};
        """

rule calculate_downsample_tr:
    """
    Compute the size to downsample each file to. If the config.yaml process parameter 'downsample'
    for the process being run is to_min_pairs_size, a group of .pairs files to compare must be
    defined. Within that group, the size of each .pairs file is computed, and the ratio of
    [smallest .pair]/[current .pair] is computed. The current .pair file will be downsampled by that
    ratio.

    Alternatively, all files can be downsampled to a constant ratio. To leave them unchanged, use '1'
    or leave the downsample parameter unspecified.
    """
    input:
        lambda wc: tree_expand("results/{proc}/{cond}/{br}/{tr}/{tr}_dedup.pairs",
                                        config['processes'][wc.proc],
                                        [('conditions', 'cond'), ('bioreps', 'br'), ('techreps', 'tr')],
                                        proc = wc.proc)
    output: "results/{proc}/.markers/calculate_downsample_tr"
    conda: "../envs/minimal_python.yaml"
    params:
        process = lambda wc: config['processes'][wc.proc],
        input_type='techreps'
    benchmark: "benchmark/downsample_tr/{proc}.tsv"
    log: "log/downsample_tr/{proc}.log"
    script: "../scripts/calculate_downsample.py"

rule downsample_tr:
    """
    Downsample all files to sizes computed in calculate_downsample_tr
    """
    input:
        pairs = "results/{proc}/{cond}/{br}/{tr}/{tr}_dedup.pairs",
        fractions = "results/{proc}/.markers/calculate_downsample_tr"
    output: "results/{proc}/{cond}/{br}/{tr}/{tr}_downsample.pairs"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_br/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/downsample_br/{proc}/{cond}/{br}/{tr}.log"
    script: "../scripts/downsample.py"

rule downsample_stats_tr:
    """
    Compute statistics on downsampled files for inclusion in multiqc report.
    """
    input: pairs = "results/{proc}/{cond}/{br}/{tr}/{tr}_dedup.pairs"
    output: "results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_downsample_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_stats/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/downsample_stats/{proc}/{cond}/{br}/{tr}.log"
    shell: """pairtools stats -o {output} {input}"""


######################################################################
# Called only if producing contact matrices from technical replicates
######################################################################

rule hic_tr:
    """
    Produce .hic files from downsampled .pairs files at specified resolutions.
    If restriction sites are provided, a per-fragments resolution will be included as well.

    Includes normalizations for juicer tools pre
    https://github.com/aidenlab/juicer/wiki/Pre

    Note that this will only be produced if explicitly required by rule all.
    """
    input:
        pairs = "results/{proc}/{cond}/{br}/{tr}/{tr}_downsample.pairs",
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'],
        juicer_tools = "resources/juicer_tools"
    output: "results/{proc}/{cond}/{br}/{tr}/{tr}.hic"
    params:
        digest    = lambda wc:  "-f" if config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'] else "",
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        norms     = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['norms'], ",", "-k"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda:
        "../envs/juicer_tools.yaml"
    benchmark: "benchmark/hic_tr/{proc}/{cond}/{br}/{tr}.tsv"
    log: "log/hic_tr/{proc}/{cond}/{br}/{tr}.log"
    shell:
        """java -Xmx20g -jar {input.juicer_tools} pre \
                    {params.digest} {input.re_sites} \
                    {params.bin_sizes} \
                    {params.norms} \
                    {input.pairs:q} \
                    {output:q} \
                    {params.assembly} \
                    &> {log}"""

rule mcool_tr:
    """
    Produce .cool file from downsampled .pairs files at highest resolution specified and
    a fragment-scale resolution if specified. Then produces a .mcool file at all specified resolutions.
    If restriction sites are provided, a per-fragments resolution will be included as well.

    Note that this will only be produced if explicitly required by rule all.
    """
    input:
        pairs = "results/{proc}/{cond}/{br}/{tr}/{tr}_downsample.pairs",
        chromsizes = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['chromsizes'],
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['cooltools']
    output:
        cool = "results/{proc}/{cond}/{br}/{tr}/{tr}.cool",
        mcool = "results/{proc}/{cond}/{br}/{tr}/{tr}.mcool",
        cool_digest = "results/{proc}/{cond}/{br}/{tr}/{tr}_digest.cool",
        mcool_digest = "results/{proc}/{cond}/{br}/{tr}/{tr}_digest.mcool"
    params:
        min_bin_size = lambda wc: min(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes']),
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda: "../envs/cooler.yaml"
    log: "log/mcool_tr/{proc}/{cond}/{br}/{tr}.log"
    benchmark: "benchmark/mcool_tr/{proc}/{cond}/{br}/{tr}.tsv"
    shell:
        """
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.chromsizes}:{params.min_bin_size} {input.pairs} {output.cool} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool} &>> {log};
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.re_sites} {input.pairs} {output.cool_digest} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool_digest} &>> {log}
        """
