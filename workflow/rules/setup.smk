rule download_index_reference:
    output:
        "resources/{0}_no_alts.fa".format(config['assembly'])
    params:
        url = config['assembly_urls'][config['assembly']]
    shell:
        """
        if [ ! -f {output} ]; then
            wget -O {output} {params.url}
            bwa index {output}
        fi
        """