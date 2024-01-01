import itertools as it

rule fastqc:
    params:
        output_dir = "results/reports/fastqc"
    input:
        "fastq/{replicate}_{read}.fq.gz"
    output:
        "results/reports/fastqc/{replicate}_{read}_fastqc.html"
    log:
        "logs/fastqc/{replicate}_{read}.log"
    shell:
        """fastqc -o {params.output_dir} {input} 2> {log}"""

rule multiqc:
    input:
        lambda wildcards: MERGES.format_template( \
                        "results/{replicate}/pairs", \
                        merge = [wildcards.merge])
    params:
        directory = "results/reports/multiqc",
        name = "{merge}_multiqc_report.html"
    output:
        name = "results/reports/multiqc/{merge}_multiqc_report.html"
    log:
        "logs/multiqc/{merge}.log"
    shell:
        """multiqc -o {params.directory} -n {params.name} {input} &> {log}"""

rule hicrep:
    input:
        "results/{node1}/matrix/{node1}_ds{downsample}.mcool",
        "results/{node2}/matrix/{node2}_ds{downsample}.mcool"
    output:
        scores = "results/reports/hicrep/scores/{node1}_{node2}_ds{downsample}.txt"
    log:
        "logs/hicrep/{node1}_{node2}_ds{downsample}.log"
    shell:
        """
        hicrep {input} {output.scores} --binSize 250000 --h 1 --dBPMax 6000000 2> {log};
        mkdir -p results/reports/hicrep/captions;
        echo results/reports/hicrep/scores/{wildcards.node1}_{wildcards.node2}_ds{wildcards.downsample}.txt {wildcards.downsample} {wildcards.node1} {wildcards.node2} >> results/reports/hicrep/captions/{wildcards.downsample}.txt
        """

rule hicrep_report:
    input:
        scores = lambda wildcards: expand("results/reports/hicrep/scores/{node1}_{node2}_ds{downsample}.txt", \
                            node1 = MERGES.nodes(), \
                            node2 = MERGES.nodes(), \
                            downsample = wildcards.downsample)
    output:
        "results/reports/hicrep/hicrep_ds{downsample}.png"
    log:
        "logs/hicrep_report/ds{downsample}.log"
    shell:
        """
        python workflow/scripts/hicrep_figure.py --scores {input.scores} --caption results/reports/hicrep/captions/{wildcards.downsample}.txt --output {output} &> {log}
        """
