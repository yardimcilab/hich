rule mustache:
    input:
        "results/{node1}/matrix/{node1}_ds{downsample}.hic"
    output:
        "results/{node1}/loops/{node1}_ds{downsample}_mustache.tsv"
    params:
        resolution = config["mustache"]["resolution"],
        chromosomes = config["mustache"]["chromosomes"],
        pThreshold = config["mustache"]["pThreshold"],
        sigmaZero = config["mustache"]["sigmaZero"],
        sparsityThreshold = config["mustache"]["sparsityThreshold"],
        normalization = config["mustache"]["normalization"],
    log:
        "logs/mustache/{replicate}_ds{downsample}.log"
    conda:
        "../envs/mustache.yaml"
    shell:
        """mustache {params.resolution} \
                    {params.chromosomes} \
                    {params.pThreshold} \
                    {params.sigmaZero} \
                    {params.sparsityThreshold} \
                    {params.normalization} \
                    -f {input} \
                    -o {output}"