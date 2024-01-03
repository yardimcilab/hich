import itertools as it

rule fastqc:
    params:
        output_dir = "results/reports/fastqc"
    input:
        "fastq/{replicate}_{read}.fq.gz"
    output:
        "results/reports/fastqc/{replicate}_{read}_fastqc.html"
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/{replicate}_{read}.log"
    shell:
        """fastqc -o {params.output_dir} {input} 2> {log}"""

rule multiqc:
    input:
        dir = lambda wildcards: MERGES.format_template( \
                        "results/{replicate}/pairs/", \
                        merge = [wildcards.merge]),
        parse = lambda wildcards: MERGES.format_template( \
                        "results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt", \
                        merge = [wildcards.merge]),
        dedup = lambda wildcards: MERGES.format_template( \
                        "results/{replicate}/pairs/{replicate}_ds{downsample}_sort_dedup_stats.txt", \
                        merge = [wildcards.merge], downsample = config["downsample"]),
        
    params:
        directory = "results/reports/multiqc",
        name = "{merge}_multiqc_report.html"
    output:
        name = "results/reports/multiqc/{merge}_multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc/{merge}.log"
    shell:
        """multiqc -m pairtools -f -o {params.directory} -n {params.name} {input.dir} &> {log}"""

rule hicrep:
    input:
        "results/{node1}/matrix/{node1}_ds{downsample}.mcool",
        "results/{node2}/matrix/{node2}_ds{downsample}.mcool"
    output:
        scores = "results/reports/hicrep/scores/{node1}_{node2}_ds{downsample}.txt",
    params:
        binSize = config["hicrep"]["binSize"],
        h = config["hicrep"]["h"],
        dBPMax = config["hicrep"]["dBPMax"],
        bDownSample = config["hicrep"]["bDownSample"],
        chrNames = config["hicrep"]["chrNames"],
        excludeChr = config["hicrep"]["excludeChr"]

    conda:
        "../envs/hicrep.yaml"
    log:
        "logs/hicrep/{node1}_{node2}_ds{downsample}.log"
    shell:
        """
        hicrep {input} {output.scores} {params.binSize} {params.h} {params.dBPMax} {params.bDownSample} {params.chrNames} {params.excludeChr} 2> {log};
        mkdir -p results/reports/hicrep/captions;
        echo results/reports/hicrep/scores/{wildcards.node1}_{wildcards.node2}_ds{wildcards.downsample}.txt {wildcards.downsample} {wildcards.node1} {wildcards.node2} >> results/reports/hicrep/captions/{wildcards.downsample}.txt
        """

rule hicrep_heatmap:
    input:
        scores = lambda wildcards: expand("results/reports/hicrep/scores/{node1}_{node2}_ds{downsample}.txt", \
                            node1 = list(MERGES), \
                            node2 = list(MERGES), \
                            downsample = wildcards.downsample)
    output:
        "results/reports/hicrep/hicrep_ds{downsample}.png"
    conda:
        "../envs/hicrep_heatmap.yaml"
    log:
        "logs/hicrep_heatmap/ds{downsample}.log"
    shell:
        """
        python workflow/scripts/hicrep_heatmap.py --scores {input.scores} --caption results/reports/hicrep/captions/{wildcards.downsample}.txt --output {output} &> {log}
        """
