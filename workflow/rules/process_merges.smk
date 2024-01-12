rule pairtools_merge:
    input:
        lambda wildcards: MERGES.format_template( \
            "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_select.pairs", \
            downsample = [wildcards.downsample],
            merge = [wildcards.merge])
    output:
        "results/{merge}/pairs/{merge}_ds{downsample}.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_merge/{merge}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/pairtools_merge/{merge}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        """pairtools merge -o {output:q} {input:q} &> {log}"""

rule pairtools_sort_merge:
    input:
        "results/{merge}/pairs/{merge}_ds{downsample}.pairs"
    output:
        "results/{merge}/pairs/{merge}_ds{downsample}_sort.pairs"
    conda:
        "../envs/pairtools.yaml"
    log:
        "logs/pairtools_sort_merge/{merge}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/pairtools_sort_merge/{merge}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        """pairtools sort -o {output:q} {input:q} &> {log}"""

rule hic_merge:
    input:
        "results/{merge}/pairs/{merge}_ds{downsample}_sort.pairs"
    output:
        "results/{merge}/matrix/{merge}_ds{downsample}.hic"
    params:
        resolutions = ','.join([str(resolution) for resolution in config['resolutions']]),
        genomeID = config['genomeID'],
        juicer_tools_jar = config['juicer_tools_jar'],
        restriction_sites = config['restriction_site_file']['juicer_tools_pre'],
        norms = ','.join([str(norm) for norm in config['normalizations']])
    conda:
        "../envs/juicer_tools.yaml"
    log:
        "logs/hic_merge/{merge}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/hic_merge/{merge}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre    -f {params.restriction_sites} \
                                                                -r {params.resolutions:q} \
                                                                -k {params.norms:q} \
                                                                {input:q} \
                                                                {output:q} \
                                                                {params.genomeID} \
                                                                &> {log}"""

rule mcool_merge:
    input:
        #"results/{merge}/matrix/{merge}_ds{downsample}.hic"
        binned = lambda wildcards: MERGES.format_template( \
            "results/{replicate}/matrix/{replicate}_ds{downsample}.cool", \
            downsample = [wildcards.downsample],
            merge = [wildcards.merge]),
        digest = lambda wildcards: MERGES.format_template( \
            "results/{replicate}/matrix/{replicate}_ds{downsample}.cool", \
            downsample = [wildcards.downsample],
            merge = [wildcards.merge])
    output:
        cool = "results/{merge}/matrix/{merge}_ds{downsample}.cool",
        mcool = "results/{merge}/matrix/{merge}_ds{downsample}.mcool",
        cool_digest = "results/{merge}/matrix/{merge}_ds{downsample}_digest.cool",
        mcool_digest = "results/{merge}/matrix/{merge}_ds{downsample}_digest.mcool"
    params:
        assembly = config["assembly_header"],
        binsizes = ','.join([str(r) for r in config["resolutions"]]),
        chromsizes = config["chromsizes"]
    conda:
        #"../envs/hic2cool.yaml"
        "../envs/cooler.yaml"
    log:
        "logs/mcool_merge/{merge}_ds{downsample}.log"
    benchmark:
        repeat("benchmarks/mcool_merge/{merge}_ds{downsample}.tsv", config["benchmark_repeat"])
    shell:
        #"""hic2cool convert {input:q} {output:q} -p 24 &> {log}"""
        # Still need to do cooler digest
        """
        cooler merge {output.cool} {input.binned};
        cooler zoomify -r {params.binsizes} --balance --balance-args '--max-iters 1000' {output.cool} &> {log};
        cooler merge {output.cool_digest} {input.digest};
        cooler zoomify -r {params.binsizes} --balance --balance-args '--max-iters 1000' {output.cool_digest} &> {log}
        """