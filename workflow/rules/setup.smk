rule download_index_reference:
    params:
        url = config['assembly_urls'][config['assembly_header']]
    output:
        config['genome_prefix'] + ".fa"
    log:
        wget = "logs/download_index_reference/wget.log",
        bwa_index = "logs/download_index_reference/bwa_index.log",
    shell:
        """
        if [ ! -f {output} ]; then
            wget -O {output} {params.url} 2> {log.wget}
            bwa index {output} 2> {log.bwa_index}
        fi
        """