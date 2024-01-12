include: "merge_plan.smk"

rule align:
    input:
        r1 = "fastq/{replicate}_1.fq.gz",
        r2 = "fastq/{replicate}_2.fq.gz"
    output:
        "results/{replicate}/sambam/{replicate}.bam"
    params:
        genome_prefix = config['genome_prefix'],
        min_mapq = config['bwa_mem_min_mapq']
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        align = "logs/align/{replicate}.log",
        to_bam = "logs/to_bam/{replicate}.log"
    benchmark:
        repeat("benchmarks/align/{replicate}.tsv", config["benchmark_repeat"])
    resources:
        mem_gb=10
    shell:
        """bwa mem -SP5M {params.min_mapq} -P {params.genome_prefix:q} {input.r1:q} {input.r2:q} 2> {log.align} | samtools view -b -o {output:q} 2> {log.to_bam}"""

rule name_sort:
    input:
        "results/{replicate}/sambam/{replicate}.bam"
    output:
        "results/{replicate}/sambam/{replicate}.name_sort.bam"
    conda:
        "../envs/bwa_samtools.yaml"
    log:
        "logs/name_sort/{replicate}.log"
    benchmark:
        repeat("benchmarks/name_sort/{replicate}.tsv", config["benchmark_repeat"])
    shell:
        """samtools sort -n -o {output:q} {input:q} 2> {log}"""

rule pairtools_parse:
    input:
        name_sort = "results/{replicate}/sambam/{replicate}.name_sort.bam"
    output:
        pairtools_parse = "results/{replicate}/pairs/{replicate}.pairs",
        pairtools_parse_stats = "results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt",
        pairtools_directory = directory("results/{replicate}/pairs")
    params:
        assembly_header = config['assembly_header'],
        min_mapq = config['pairtools_parse_min_mapq'],
        chromsizes = config['chromsizes']
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_parse/{replicate}.log"
    threads: 24
    benchmark:
        repeat("benchmarks/pairtools_parse/{replicate}.tsv", config["benchmark_repeat"])
    shell:
        """
        pairtools parse {input:q} \
            {params.min_mapq} \
            -o {output.pairtools_parse:q} \
            --output-stats {output.pairtools_parse_stats:q} \
            -c {params.chromsizes:q} \
            --assembly {params.assembly_header} \
            --nproc-in {threads} \
            --nproc-out {threads} \
            &> {log}
        """

rule wait_parse_all:
    input:
        [MERGES.format_template("results/{replicate}/pairs/{replicate}.pairs")]
    output:
        "wait_parse_all"
    log:
        "logs/wait_parse_all/wait_parse_all.log"
    benchmark:
        repeat("benchmarks/wait_parse_all/wait_parse_all.tsv", config["benchmark_repeat"])
    shell:
        """touch {output:q} &> {log}"""

rule pairtools_downsample:
    input:
        parse = "results/{replicate}/pairs/{replicate}.pairs",
        marker = "wait_parse_all"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_downsample/{replicate}_ds{downsample}.log"
    params:
        replicates = " ".join(sorted(list(MERGES.source_nodes())))
    benchmark:
        repeat("benchmarks/pairtools_downsample/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
    script:
        "../scripts/downsample.py"

rule pairtools_sort:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}.pairs"
    output:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_sort/{replicate}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/pairtools_sort/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
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
    benchmark:
        repeat("benchmarks/pairtools_deduplicate/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
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
    benchmark:
        repeat("benchmarks/pairtools_select/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
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
        restriction_sites = config['restriction_site_file']['juicer_tools_pre'],
        norms = ','.join([str(norm) for norm in config['normalizations']])
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs"
    output:
        "results/{replicate}/matrix/{replicate}_ds{downsample}.hic"
    conda:
        "../envs/juicer_tools.yaml"
    log:
        "logs/hic/{replicate}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/hic/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre \
                    -f {params.restriction_sites} \
                    -r {params.resolutions:q} \
                    -k {params.norms:q} \
                    {input:q} \
                    {output:q} \
                    {params.genomeID} \
                    &> {log}"""

rule mcool:
    input:
        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup.pairs"
    output:
        cool = "results/{replicate}/matrix/{replicate}_ds{downsample}.cool",
        mcool = "results/{replicate}/matrix/{replicate}_ds{downsample}.mcool",
        cool_digest = "results/{replicate}/matrix/{replicate}_ds{downsample}_digest.cool",
        mcool_digest = "results/{replicate}/matrix/{replicate}_ds{downsample}_digest.mcool"
    params:
        assembly = config["assembly_header"],
        min_binsize = min(config["resolutions"]),
        binsizes = ','.join([str(r) for r in config["resolutions"]]),
        chromsizes = config["chromsizes"],
        restriction_sites = config['restriction_site_file']['bed_format']
    conda:
        "../envs/cooler.yaml"
    log:
        "logs/mcool/{replicate}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/mcool/{replicate}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        """
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {params.chromsizes}:{params.min_binsize} {input} {output.cool} &> {log};
        cooler zoomify -r {params.binsizes} --balance --balance-args '--max-iters 1000' {output.cool} &> {log};
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {params.restriction_sites} {input} {output.cool_digest} &> {log};
        cooler zoomify -r {params.binsizes} --balance --balance-args '--max-iters 1000' {output.cool_digest} &> {log}
        """
