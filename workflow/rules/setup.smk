rule get_noalts_reference:
    output:
        "resources/{assembly}_no_alts.fa"
    params:
        urls = {"hg38":"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"}
    shell:
        """
        if [ ! -f {output} ]; then
            wget -O {output} {params.urls[config['assembly']]}
        fi
        """