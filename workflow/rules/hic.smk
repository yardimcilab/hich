include: "setup.smk"

rule align_hic:
    params:
        ref = config['assembly']
    input:
        r1 = "fastq/{replicate}_1.fq.gz",
        r2 = "fastq/{replicate}_2.fq.gz"
    output:
        "results/{experiment}/sambam/{replicate}.bam"
    shell:
        """bwa mem -SP5M -t 24 {params.ref:q} {input.r1:q} {input.r2:q} | samtools view -b -o {output:q}"""

rule name_sort:
    input:
        "results/{experiment}/sambam/{replicate}.bam"
    output:
        "results/{experiment}/sambam/{replicate}.name_sort.bam"
    shell:
        """samtools sort -n -o {output:q} {input:q}"""

rule pairtools_parse:
    params:
        assembly = config['assembly'],
        min_mapq = config['min_mapq'],
        chromsizes = config['chromsizes']
    input:
        name_sort = "results/{experiment}/sambam/{replicate}.name_sort.bam"
    output:
        pairtools_parse = "results/{experiment}/{replicate}/pairs/{replicate}.pairs",
        pairtools_parse_stats = "results/{experiment}/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt"
    
    shell:
        """
        pairtools parse {input:q} \
            -o {output.pairtools_parse:q} \
            --output-stats {output.pairtools_parse_stats:q} \
            -c {params.chromsizes:q} \
            --assembly {params.assembly} \
            --min-mapq {params.min_mapq} \
            --nproc-in 24 \
            --nproc-out 24
        """

rule wait_parse_all:
    input:
        [filename for filename, _, _ in EXPERIMENTS.FormatTemplate("results/{experiment}/{replicate}/pairs/{replicate}.pairs")]
    output:
        "wait_parse_all"
    shell:
        """touch {output:q}"""

rule pairtools_downsample:
    input:
        parse = "results/{experiment}/{replicate}/pairs/{replicate}.pairs",
        marker = "wait_parse_all"
    output:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    run:
        downsample = wildcards.downsample
        if downsample == "min":
            total_mapped = compute_replicate_total_mapped("results/{experiment}/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt")
            downsample = min_downsample(total_mapped, wildcards)

        shell(f"""pairtools sample --output \"{output}\" --seed 0 {downsample} \"{input.parse}\"""")

rule pairtools_sort:
    input:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    output:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    shell:
        """pairtools sort -o {output:q} {input:q}"""

rule pairtools_deduplicate:
    input:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    output:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs"
    shell:
        """pairtools dedup --mark-dups -o {output:q} {input:q}"""

rule pairtools_select:
    input:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs"
    output:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs"
    shell:
        """pairtools select '(pair_type=="UU") or (pair_type=="RU") or (pair_type=="UR")' -o {output:q} {input:q}"""

rule hic:
    params:
        resolutions = config['resolutions'],
        assembly = config['assembly'],
        juicer_tools_jar = config['juicer_tools_jar']
    input:
        "results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs"
    output:
        "results/{experiment}/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre -r {params.resolutions:q} {input:q} {output:q} {params.assembly}"""

rule mcool:
    input:
        "results/{experiment}/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    output:
        "results/{experiment}/{replicate}/matrix/{replicate}_ds{downsample}.mcool"
    shell:
        """hic2cool convert {input:q} {output:q} -p 24"""

rule pairtools_merge:
    input:
        lambda wildcards: [filename for filename, _, _ in EXPERIMENTS.FormatTemplate \
            ("results/{experiment}/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs", \
            downsample = wildcards.downsample)]
    output:
        "results/{experiment}/pairs/{experiment}_ds{downsample}_sort_dedup_select_merge.pairs"
    shell:
        """pairtools merge -o {output:q} {input:q}"""

rule pairtools_sort_merge:
    input:
        "results/{experiment}/pairs/{experiment}_ds{downsample}_sort_dedup_select_merge.pairs"
    output:
        "results/{experiment}/pairs/{experiment}_ds{downsample}_sort_dedup_select_merge.pairs"
    shell:
        """pairtools sort -o {output:q} {input:q}"""

rule hic_merge:
    params:
        resolutions = config['resolutions'],
        assembly = config['assembly'],
        juicer_tools_jar = config['juicer_tools_jar']
    input:
        "results/{experiment}/pairs/{experiment}_ds{downsample}_sort_dedup_select_merge.pairs"
    output:
        "results/{experiment}/matrix/{experiment}_ds{downsample}.hic"
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre -r {params.resolutions:q} {input:q} {output:q} {params.assembly}"""

rule mcool_merge:
    input:
        "results/{experiment}/matrix/{experiment}_ds{downsample}.hic"
    output:
        "results/{experiment}/matrix/{experiment}_ds{downsample}.mcool"
    shell:
        """hic2cool convert {input:q} {output:q} -p 24"""