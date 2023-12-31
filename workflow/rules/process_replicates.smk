include: "merge_plan.smk"

rule align:
    input:
        r1 = "fastq/{replicate}_1.fq.gz",
        r2 = "fastq/{replicate}_2.fq.gz"
    output:
        "results/{replicate}/sambam/{replicate}.bam"
    params:
        genome_prefix = config['genome_prefix']
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        align = "logs/align/{replicate}.log",
        to_bam = "logs/to_bam/{replicate}.log"
    shell:
        """bwa mem -SP5M -t 24 -P {params.genome_prefix:q} {input.r1:q} {input.r2:q} 2> {log.align} | samtools view -b -o {output:q} 2> {log.to_bam}"""

rule name_sort:
    input:
        "results/{replicate}/sambam/{replicate}.bam"
    output:
        "results/{replicate}/sambam/{replicate}.name_sort.bam"
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        "logs/name_sort/{replicate}.log"
    shell:
        """samtools sort -n -o {output:q} {input:q} 2> {log}"""

rule pairtools_parse:
    input:
        name_sort = "results/{replicate}/sambam/{replicate}.name_sort.bam"
    output:
        pairtools_parse = "results/{replicate}/pairs/{replicate}.pairs",
        pairtools_parse_stats = "results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt"
    params:
        assembly_header = config['assembly_header'],
        min_mapq = config['min_mapq'],
        chromsizes = config['chromsizes']
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_parse/{replicate}.log"
    shell:
        """
        pairtools parse {input:q} \
            -o {output.pairtools_parse:q} \
            --output-stats {output.pairtools_parse_stats:q} \
            -c {params.chromsizes:q} \
            --assembly {params.assembly_header} \
            --nproc-in 24 \
            --nproc-out 24 \
            &> {log}
        """

rule wait_parse_all:
    input:
        [MERGES.format_template("results/{replicate}/pairs/{replicate}.pairs")]
    output:
        "wait_parse_all"
    log:
        "logs/wait_parse_all/wait_parse_all.log"
    shell:
        """touch {output:q} &> {log}"""

rule pairtools_downsample:
    input:
        parse = "results/{replicate}/pairs/{replicate}.pairs",
        marker = "wait_parse_all"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    log:
        "logs/pairtools_downsample/{replicate}_ds{downsample}.log"
    run:
        downsample = wildcards.downsample
        if downsample == "min":
            total_mapped = compute_replicate_total_mapped("results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt", MERGES)
            downsample = min_downsample(total_mapped, wildcards)

        shell(f"""pairtools sample \
                            --output \"{output}\" \
                            --seed 0 \
                            {downsample} \
                            \"{input.parse}\" \
                            &> {log}""")

rule pairtools_sort:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_sort/{replicate}_ds{downsample}.log"
    shell:
        """pairtools sort \
                    -o {output:q} \
                    {input:q} \
                    &> {log}"""

rule pairtools_deduplicate:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    output:
        pairtools_deduplicate = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs",
        pairtools_duplicates = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_removed_duplicates.pairs",
        pairtools_deduplicate_stats = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_stats.txt"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_deduplicate/{replicate}_ds{downsample}.log"
    shell:
        """pairtools dedup \
                    --send-header-to both \
                    --mark-dups \
                    --output-stats {output.pairtools_deduplicate_stats:q} \
                    --mark-dups \
                    --output-dups {output.pairtools_duplicates:q} \
                    -o {output.pairtools_deduplicate:q} \
                    {input:q} \
                    &> {log}"""

rule pairtools_select:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs"
    output:
        pairtools_select = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs",
        pairtools_rest = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_rest.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_select/{replicate}_ds{downsample}.log"
    shell:
        """pairtools select '(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")' \
                    --output-rest {output.pairtools_rest} \
                    -o {output.pairtools_select:q} \
                    {input:q} \
                    &> {log}"""

rule hic:
    params:
        resolutions = ','.join([str(resolution) for resolution in config['resolutions']]),
        genomeID = config['genomeID'],
        juicer_tools_jar = config['juicer_tools_jar'],
        min_mapq = config['min_mapq'],
        restriction_sites = config['juicer_tools_pre_restriction_site_file'],
        norms = ','.join([str(norm) for norm in config['normalizations']])
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs"
    output:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    conda:
        "../envs/juicer_tools.yaml"
    log:
        "logs/hic/{replicate}_ds{downsample}.log"
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre \
                    -q {params.min_mapq} \
                    -f {params.restriction_sites} \
                    -r {params.resolutions:q} \
                    -k {params.norms:q} \
                    {input:q} \
                    {output:q} \
                    {params.genomeID} \
                    &> {log}"""

rule mcool:
    input:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    output:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.mcool"
    conda:
        "../envs/hic2cool.yaml"
    log:
        "logs/mcool/{replicate}_ds{downsample}.log"
    shell:
        """hic2cool convert {input:q} {output:q} -p 24 &> {log}"""