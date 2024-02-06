rule merge_br:
    """
    Merge biological replicate .pairs files into a single condition .pairs and
    produce output statistics for the conditions multiqc report
    """
    input:
        pairs = lambda wc: expand("results/{proc}/{cond}/{br}/{br}_downsample.pairs", \
                        cond = wc.cond,
                        proc = wc.proc,
                        br = config['processes'][wc.proc]['conditions'][wc.cond]['bioreps']),
        stats = lambda wc: expand("results/{proc}/reports/multiqc/bioreps/{cond}_{br}_downsample_stats.txt",
                        proc = wc.proc,
                        cond = wc.cond,
                        br = config['processes'][wc.proc]['conditions'][wc.cond]['bioreps'])
    output:
        pairs = "results/{proc}/{cond}/{cond}.pairs",
        stats = "results/{proc}/reports/multiqc/conditions/{cond}_merge_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/merge/{proc}/{cond}.log"
    log: "log/merge/{proc}/{cond}.log"
    shell: """pairtools merge --keep-first-header -o {output.pairs:q} {input.pairs:q} &>> {log}; pairtools stats --merge -o {output.stats:q} {input.stats:q};"""

rule calculate_downsample_cd:
    """
    Compute the size to downsample each condition to. If the config.yaml process parameter 'downsample'
    for the process being run is to_min_pairs_size, a group of .pairs files to compare must be
    defined. Within that group, the size of each .pairs file is computed, and the ratio of
    [smallest .pair]/[current .pair] is computed. The current .pair file will be downsampled by that
    ratio.

    Alternatively, all files can be downsampled to a constant ratio. To leave them unchanged, use '1'
    or leave the downsample parameter unspecified.
    """
    input:
        lambda wc: tree_expand("results/{proc}/{cond}/{cond}.pairs",
                                        config['processes'][wc.proc],
                                        [('conditions', 'cond'), ('bioreps', 'br')],
                                        proc = wc.proc)
    output: "results/{proc}/.markers/calculate_downsample_cd"
    conda: "../envs/minimal_python.yaml"
    benchmark: "benchmark/calculate_downsample_br/{proc}.tsv"
    log: "log/calculate_downsample_br/{proc}.log"
    params:
        process = lambda wc: config['processes'][wc.proc],
        input_type='conditions'
    script: "../scripts/calculate_downsample.py"

rule downsample_cd:
    """
    Downsample all conditions to sizes computed in calculate_downsample_cd
    """
    input:
        pairs = "results/{proc}/{cond}/{cond}.pairs",
        fractions = "results/{proc}/.markers/calculate_downsample_cd"
    output: "results/{proc}/{cond}/{cond}_downsample.pairs"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_br/{proc}/{cond}.tsv"
    log: "log/downsample_br/{proc}/{cond}.log"
    script: "../scripts/downsample.py"

rule downsample_stats_cd:
    """
    Compute statistics on downsampled conditions for inclusion in multiqc report.
    """
    input:
        pairs = "results/{proc}/{cond}/{cond}_downsample.pairs"
    output: "results/{proc}/reports/multiqc/conditions/{cond}_downsample_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_stats/{proc}/{cond}.tsv"
    log: "log/downsample_stats/{proc}/{cond}.log"
    shell: """pairtools stats -o {output} {input}"""

rule hic_cd:
    """
    Produce .hic files from downsampled single condition .pairs files at specified resolutions.
    If restriction sites are provided, a per-fragments resolution will be included as well.

    Includes normalizations for juicer tools pre
    https://github.com/aidenlab/juicer/wiki/Pre

    Note that this will only be produced if explicitly required by rule all.
    """
    input:
        pairs = "results/{proc}/{cond}/{cond}_downsample.pairs",
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'],
        juicer_tools = "resources/juicer_tools"
    output: "results/{proc}/{cond}/{cond}.hic"
    params:
        digest    = lambda wc:  "-f" if config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'] else "",
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        norms     = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['norms'], ",", "-k"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda:
        "../envs/juicer_tools.yaml"
    benchmark: "benchmark/hic/{proc}/{cond}/{cond}.tsv"
    log: "log/hic/{proc}/{cond}/{cond}.log"
    shell:
        """java -Xmx20g -jar {input.juicer_tools} pre \
                    {params.digest} {input.re_sites} \
                    {params.bin_sizes} \
                    {params.norms} \
                    {input.pairs:q} \
                    {output:q} \
                    {params.assembly} \
                    &> {log}"""

rule mcool_cd:
    """
    Produce .cool file from downsampled single condition .pairs files at highest resolution specified and
    a fragment-scale resolution if specified. Then produces a .mcool file at all specified resolutions.
    If restriction sites are provided, a per-fragments resolution will be included as well.

    Note that this will only be produced if explicitly required by rule all.
    """
    input:
        pairs = "results/{proc}/{cond}/{cond}_downsample.pairs",
        chromsizes = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['chromsizes'],
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['cooltools']
    output:
        cool = "results/{proc}/{cond}/{cond}.cool",
        mcool = "results/{proc}/{cond}/{cond}.mcool",
        cool_digest = "results/{proc}/{cond}/{cond}_digest.cool",
        mcool_digest = "results/{proc}/{cond}/{cond}_digest.mcool"
    params:
        min_bin_size = lambda wc: min(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes']),
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda: "../envs/cooler.yaml"
    log: "log/mcool/{proc}/{cond}/{cond}.log"
    benchmark: "benchmark/mcool/{proc}/{cond}/{cond}.tsv"
    shell:
        """
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.chromsizes}:{params.min_bin_size} {input.pairs} {output.cool} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool} &>> {log};
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.re_sites} {input.pairs} {output.cool_digest} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool_digest} &>> {log}
        """