# Rewire to take from downsample_tr
rule merge_tr:
    input:
        pairs = lambda wc: expand("results/{proc}/{cond}/{br}/{tr}/{tr}_downsample.pairs",
                        proc = wc.proc,
                        cond = wc.cond,
                        br = wc.br,
                        tr = config['processes'][wc.proc]['conditions'][wc.cond]['bioreps'][wc.br]['techreps']),
        stats = lambda wc: expand("results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_downsample_stats.txt",
                        proc = wc.proc,
                        cond = wc.cond,
                        br = wc.br,
                        tr = config['processes'][wc.proc]['conditions'][wc.cond]['bioreps'][wc.br]['techreps'])
    output:
        pairs = "results/{proc}/{cond}/{br}/{br}.pairs",
        stats = "results/{proc}/reports/multiqc/bioreps/{cond}_{br}_merge_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/merge_br/{proc}/{cond}/{br}.log"
    log: "log/merge_br/{proc}/{cond}/{br}.log"
    shell: """pairtools merge -o {output.pairs:q} {input.pairs:q} &>> {log}; pairtools stats --merge -o {output.stats:q} {input.stats:q};"""

rule dedup_br:
    input: "results/{proc}/{cond}/{br}/{br}.pairs"
    output:
        dedup = "results/{proc}/{cond}/{br}/{br}_dedup.pairs",
        dups = "results/{proc}/{cond}/{br}/{br}_duplicates.pairs",
        stats = "results/{proc}/reports/multiqc/bioreps/{cond}_{br}_dedup_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/dedup_br/{proc}/{cond}/{br}.tsv"
    log: "log/dedup_br/{proc}/{cond}/{br}.log"
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

rule calculate_downsample_br:
    input:
        lambda wc: tree_expand("results/{proc}/{cond}/{br}/{br}_dedup.pairs",
                                        config['processes'][wc.proc],
                                        [('conditions', 'cond'), ('bioreps', 'br')],
                                        proc = wc.proc)
    output: "results/{proc}/.markers/calculate_downsample_br"
    benchmark: "benchmark/calculate_downsample_br/{proc}.tsv"
    log: "log/calculate_downsample_br/{proc}.log"
    params:
        process = lambda wc: config['processes'][wc.proc],
        input_type='bioreps'

    script: "../scripts/calculate_downsample.py"

rule downsample_br:
    input:
        pairs = "results/{proc}/{cond}/{br}/{br}_dedup.pairs",
        fractions = "results/{proc}/.markers/calculate_downsample_br"
    output: "results/{proc}/{cond}/{br}/{br}_downsample.pairs"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_br/{proc}/{cond}/{br}.tsv"
    log: "log/downsample_br/{proc}/{cond}/{br}.log"
    script: "../scripts/downsample.py"

rule downsample_stats_br:
    input:
        pairs = "results/{proc}/{cond}/{br}/{br}_downsample.pairs"
    output: "results/{proc}/reports/multiqc/bioreps/{cond}_{br}_downsample_stats.txt"
    container: "docker://bskubi/pairtools:1.0.3"
    benchmark: "benchmark/downsample_stats/{proc}/{cond}/{br}.tsv"
    log: "log/downsample_stats/{proc}/{cond}/{br}.log"
    shell: """pairtools stats -o {output} {input}"""

rule hic_br:
    input:
        pairs = "results/{proc}/{cond}/{br}/{br}_downsample.pairs",
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'],
        juicer_tools = "resources/juicer_tools"
    output: "results/{proc}/{cond}/{br}/{br}.hic"
    params:
        digest    = lambda wc:  "-f" if config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['juicer_tools'] else "",
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        norms     = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['norms'], ",", "-k"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda:
        "../envs/juicer_tools.yaml"
    benchmark: "benchmark/hic_br/{proc}/{cond}/{br}.tsv"
    log: "log/hic_br/{proc}/{cond}/{br}.log"
    shell:
        """java -Xmx20g -jar {input.juicer_tools} pre \
                    {params.digest} {input.re_sites} \
                    {params.bin_sizes} \
                    {params.norms} \
                    {input.pairs:q} \
                    {output:q} \
                    {params.assembly} \
                    &> {log}"""

rule mcool_br:
    input:
        pairs = "results/{proc}/{cond}/{br}/{br}_downsample.pairs",
        chromsizes = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['chromsizes'],
        re_sites = lambda wc: config['processes'][wc.proc]['conditions'][wc.cond]['genome']['restriction_sites']['cooltools']
    output:
        cool = "results/{proc}/{cond}/{br}/{br}.cool",
        mcool = "results/{proc}/{cond}/{br}/{br}.mcool",
        cool_digest = "results/{proc}/{cond}/{br}/{br}_digest.cool",
        mcool_digest = "results/{proc}/{cond}/{br}/{br}_digest.mcool"
    params:
        min_bin_size = lambda wc: min(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes']),
        bin_sizes = lambda wc:  fmtl(config['processes'][wc.proc]['parameters']['matrix']['bin_sizes'], ",", "-r"),
        assembly  = lambda wc:  config['processes'][wc.proc]['conditions'][wc.cond]['genome']['assembly']
    conda: "../envs/cooler.yaml"
    log: "log/mcool_br/{proc}/{cond}/{br}.log"
    benchmark: "benchmark/mcool_br/{proc}/{cond}/{br}.tsv"
    shell:
        """
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.chromsizes}:{params.min_bin_size} {input.pairs} {output.cool} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool} &>> {log};
        cooler cload pairs --assembly {params.assembly} -c1 2 -p1 3 -c2 4 -p2 5 {input.re_sites} {input.pairs} {output.cool_digest} &>> {log};
        cooler zoomify {params.bin_sizes} --balance --balance-args '--max-iters 1000' {output.cool_digest} &>> {log}
        """