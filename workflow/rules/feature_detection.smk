rule mustache:
    input:
        "results/{node1}/matrix/{node1}_ds{downsample}.hic"
    output:
        "results/reports/mustache/{node1}/loops/{node1}_ds{downsample}_{mustache_suffix}.tsv"
    params:
        resolution = config["mustache"]["resolution"],
        distance = config["mustache"]["distance"],
        chromosomes = config["mustache"]["chromosomes"],
        pThreshold = config["mustache"]["pThreshold"],
        sigmaZero = config["mustache"]["sigmaZero"],
        sparsityThreshold = config["mustache"]["sparsityThreshold"],
        normalization = config["mustache"]["normalization"]
    log:
        "logs/mustache/{node1}_ds{downsample}_{mustache_suffix}.log"
    conda:
        "../envs/mustache.yaml"
    shell:
        """mustache {params.resolution} \
                    {params.distance} \
                    {params.chromosomes} \
                    {params.pThreshold} \
                    {params.sigmaZero} \
                    {params.sparsityThreshold} \
                    {params.normalization} \
                    -f {input} \
                    -o {output} &> {log}"""

rule diff_mustache:
    input:
        f1 = "results/{node1}/matrix/{node1}_ds{downsample}.hic",
        f2 = "results/{node2}/matrix/{node2}_ds{downsample}.hic"
    output:
        diffloop1 = "results/reports/diff_mustache/{node1}_{node2}_ds{downsample}_{diff_mustache_suffix}.tsv.diffloop1",
        diffloop2 = "results/reports/diff_mustache/{node1}_{node2}_ds{downsample}_{diff_mustache_suffix}.tsv.diffloop2"
    params:
        resolution = config["diff_mustache"]["resolution"],
        distance = config["diff_mustache"]["distance"],
        chromosomes = config["diff_mustache"]["chromosomes"],
        pThreshold = config["diff_mustache"]["pThreshold"],
        sigmaZero = config["diff_mustache"]["sigmaZero"],
        sparsityThreshold = config["diff_mustache"]["sparsityThreshold"],
        normalization = config["diff_mustache"]["normalization"],
        output_base = "results/reports/diff_mustache/{node1}_{node2}_ds{downsample}_{diff_mustache_suffix}.tsv",
    log:
        "logs/diff_mustache/{node1}_{node2}_ds{downsample}_{diff_mustache_suffix}.log"
    conda:
        "../envs/mustache.yaml"
    shell:
        """python -m mustache.diff_mustache   {params.resolution} \
                                                                            {params.distance} \
                                                                            {params.chromosomes} \
                                                                            {params.pThreshold} \
                                                                            {params.sigmaZero} \
                                                                            {params.sparsityThreshold} \
                                                                            {params.normalization} \
                                                                            -f1 {input.f1} \
                                                                            -f2 {input.f2} \
                                                                            -o {params.output_base} &> {log};
        rm {params.output_base}.loop1 {params.output_base}.loop2
        """
