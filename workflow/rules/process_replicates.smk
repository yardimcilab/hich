include: "instantiate_merges.smk"

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
    shell:
        """bwa mem -SP5M -t 24 -P {params.genome_prefix:q} {input.r1:q} {input.r2:q} | samtools view -b -o {output:q}"""

rule name_sort:
    input:
        "results/{replicate}/sambam/{replicate}.bam"
    output:
        "results/{replicate}/sambam/{replicate}.name_sort.bam"
    conda:
        "../envs/bwa_samtools.yaml"
    shell:
        """samtools sort -n -o {output:q} {input:q}"""

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
    shell:
        """
        pairtools parse {input:q} \
            -o {output.pairtools_parse:q} \
            --output-stats {output.pairtools_parse_stats:q} \
            -c {params.chromsizes:q} \
            --assembly {params.assembly_header} \
            --min-mapq {params.min_mapq} \
            --nproc-in 24 \
            --nproc-out 24
        """

rule wait_parse_all:
    input:
        [MERGES.format_template("results/{replicate}/pairs/{replicate}.pairs")]
    output:
        "wait_parse_all"
    shell:
        """touch {output:q}"""

rule pairtools_downsample:
    input:
        parse = "results/{replicate}/pairs/{replicate}.pairs",
        marker = "wait_parse_all"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    run:
        downsample = wildcards.downsample
        if downsample == "min":
            total_mapped = compute_replicate_total_mapped("results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt", MERGES)
            downsample = min_downsample(total_mapped, wildcards)

        shell(f"""pairtools sample \
                            --output \"{output}\" \
                            --seed 0 \
                            {downsample} \
                            \"{input.parse}\"""")

rule pairtools_sort:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    conda:
        "../envs/pairtools.yaml"
    shell:
        """pairtools sort \
                    -o {output:q} \
                    {input:q}"""

rule pairtools_deduplicate:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    output:
        pairtools_deduplicate = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs",
        pairtools_duplicates = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_removed_duplicates.pairs",
        pairtools_deduplicate_stats = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_stats.txt"
    conda:
        "../envs/pairtools.yaml"
    shell:
        """pairtools dedup \
                    --send-header-to both \
                    --mark-dups \
                    --output-stats {output.pairtools_deduplicate_stats:q} \
                    --mark-dups \
                    --output-dups {output.pairtools_duplicates:q} \
                    -o {output.pairtools_deduplicate:q} \
                    {input:q}"""

rule pairtools_select:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs"
    output:
        pairtools_select = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs",
        pairtools_rest = "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_rest.pairs"
    conda:
        "../envs/pairtools.yaml"
    shell:
        """pairtools select '(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")' \
                    --output-rest {output.pairtools_rest} \
                    -o {output.pairtools_select:q} \
                    {input:q}"""

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
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre \
                    -q {params.min_mapq} \
                    -f {params.restriction_sites} \
                    -r {params.resolutions:q} \
                    -k {params.norms:q} \
                    {input:q} \
                    {output:q} \
                    {params.genomeID}"""

rule mcool:
    input:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    output:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.mcool"
    conda:
        "../envs/hic2cool.yaml"
    shell:
        """hic2cool convert {input:q} {output:q} -p 24"""