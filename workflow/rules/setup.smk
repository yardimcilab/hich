rule download_index_reference:
    params:
        url = config['assembly_urls'][config['assembly']]
    output:
        config['genome']
    shell:
        """
        if [ ! -f {output} ]; then
            wget -O {output} {params.url}
            bwa index {output}
        fi
        """