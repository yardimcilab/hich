rule fastqc:
    input:
        r1 = config['data']['fastq'] + "{cond}_{br}_{tr}_1.fq.gz",
        r2 = config['data']['fastq'] + "{cond}_{br}_{tr}_2.fq.gz"
    output:
        r1 = "results/{proc}/reports/fastqc/{cond}_{br}_{tr}_1_fastqc.html",
        r2 = "results/{proc}/reports/fastqc/{cond}_{br}_{tr}_2_fastqc.html"
    params: 
        output_dir = "results/{proc}/reports/fastqc/"
    conda: "../envs/fastqc.yaml"
    log: "log/fastqc/{proc}/{cond}/{br}/{tr}.log"
    benchmark: "benchmark/fastqc/{proc}/{cond}/{br}/{tr}.tsv"
    shell: """fastqc -o {params.output_dir} {input.r1} &> {log}; fastqc -o {params.output_dir} {input.r2} &> {log};"""

rule multiqc_tr:
    input: 
        lambda wc:  tree_expand("results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_parse_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br'), ('techreps', 'tr')],
                                    proc = wc.proc) +
                    tree_expand("results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_dedup_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br'), ('techreps', 'tr')],
                                    proc = wc.proc) +
                    tree_expand("results/{proc}/reports/multiqc/techreps/{cond}_{br}_{tr}_downsample_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br'), ('techreps', 'tr')],
                                    proc = wc.proc)
    output: "results/{proc}/reports/multiqc/techreps/techreps_multiqc_report.html"
    params:
        directory = "results/{proc}/reports/multiqc/techreps",
        name = "techreps_multiqc_report.html"
    conda: "../envs/multiqc.yaml"
    log: "log/multiqc_tr/{proc}.log"
    benchmark: "benchmark/multiqc_tr/{proc}.tsv"
    shell: """multiqc -m pairtools -f -o {params.directory} -n {params.name} {params.directory} &> {log}"""

rule multiqc_br:
    input: 
        lambda wc:  tree_expand("results/{proc}/reports/multiqc/bioreps/{cond}_{br}_merge_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br')],
                                    proc = wc.proc) + 
                    tree_expand("results/{proc}/reports/multiqc/bioreps/{cond}_{br}_dedup_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br')],
                                    proc = wc.proc) + 
                    tree_expand("results/{proc}/reports/multiqc/bioreps/{cond}_{br}_downsample_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond'), ('bioreps', 'br')],
                                    proc = wc.proc)
    output: "results/{proc}/reports/multiqc/bioreps/bioreps_multiqc_report.html"
    params:
        directory = "results/{proc}/reports/multiqc/bioreps",
        name = "bioreps_multiqc_report.html"
    conda: "../envs/multiqc.yaml"
    log: "log/multiqc_br/{proc}.log"
    benchmark: "benchmark/multiqc_br/{proc}.tsv"
    shell: """multiqc -m pairtools -f -o {params.directory} -n {params.name} {params.directory} &> {log}"""

rule multiqc_cd:
    input: 
        lambda wc:  tree_expand("results/{proc}/reports/multiqc/conditions/{cond}_merge_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond')],
                                    proc = wc.proc) +
                    tree_expand("results/{proc}/reports/multiqc/conditions/{cond}_downsample_stats.txt",
                                    config['processes'][wc.proc],
                                    [('conditions', 'cond')],
                                    proc = wc.proc)
    output: "results/{proc}/reports/multiqc/conditions/conditions_multiqc_report.html"
    params:
        directory = "results/{proc}/reports/multiqc/conditions",
        name = "conditions_multiqc_report.html"
    conda: "../envs/multiqc.yaml"
    log: "log/multiqc_cd/{proc}.log"
    benchmark: "benchmark/multiqc_cd/{proc}.tsv"
    shell: """multiqc -m pairtools -f -o {params.directory} -n {params.name} {params.directory} &> {log}"""