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
    shell:
        """pairtools merge -o {output:q} {input:q}"""

rule pairtools_sort_merge:
    input:
        "results/{merge}/pairs/{merge}_ds{downsample}.pairs"
    output:
        "results/{merge}/pairs/{merge}_ds{downsample}_sort.pairs"
    conda:
        "../envs/pairtools.yaml"
    shell:
        """pairtools sort -o {output:q} {input:q}"""

rule hic_merge:
    params:
        resolutions = ','.join([str(resolution) for resolution in config['resolutions']]),
        genomeID = config['genomeID'],
        juicer_tools_jar = config['juicer_tools_jar'],
        min_mapq = config['min_mapq'],
        restriction_sites = config['juicer_tools_pre_restriction_site_file'],
        norms = ','.join([str(norm) for norm in config['normalizations']])
    input:
        "results/{merge}/pairs/{merge}_ds{downsample}_sort.pairs"
    output:
        "results/{merge}/matrix/{merge}_ds{downsample}.hic"
    conda:
        "../envs/juicer_tools.yaml"
    shell:
        """java -Xmx20g -jar {params.juicer_tools_jar:q} pre    -q {params.min_mapq} \
                                                                -f {params.restriction_sites} \
                                                                -r {params.resolutions:q} \
                                                                -k {params.norms:q} \
                                                                {input:q} \
                                                                {output:q} \
                                                                {params.genomeID}"""

rule mcool_merge:
    input:
        "results/{merge}/matrix/{merge}_ds{downsample}.hic"
    output:
        "results/{merge}/matrix/{merge}_ds{downsample}.mcool"
    conda:
        "../envs/hic2cool.yaml"
    shell:
        """hic2cool convert {input:q} {output:q} -p 24"""