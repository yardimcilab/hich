rule install_bwa_mem2:
    output:
        "resources/bwa-mem2/bwa-mem2"
    log:
        "logs/install_bwa_mem2"
    shell:
        """
        set -e  # Stop on error
        mkdir -p resources/bwa-mem2

        # Try downloading precompiled binaries
        if curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C resources/bwa-mem2 --strip-components=1; then
            echo "Downloaded and extracted precompiled binaries successfully."
        else
            # If download fails, compile from source
            echo "Downloading precompiled binaries failed. Compiling from source."
            cd resources/bwa-mem2
            git clone --recursive https://github.com/bwa-mem2/bwa-mem2 .
            # Alternatively, for a non-recursive clone
            # git clone https://github.com/bwa-mem2/bwa-mem2 .
            # git submodule init
            # git submodule update

            make
        fi
        """

rule download_reference:
    params:
        url = config['assembly_urls'][config['assembly_header']]
    output:
        config['genome_prefix']
    log:
        "logs/download_reference/wget.log"
    shell:
        """
        wget -O {output} {params.url} 2> {log}
        """

rule index_reference:
    input:
        config['genome_prefix']
    output:
        fa = config['genome_prefix'] + ".0123",
        amb = config['genome_prefix'] + ".amb",
        ann = config['genome_prefix'] + ".ann",
        bwt = config['genome_prefix'] + ".bwt.2bit.64",
        pac = config['genome_prefix'] + ".pac"
    log:
        "logs/index_reference/index_reference.log"
    shell:
        """
        resources/bwa-mem2/bwa-mem2 index {input} &> {log}
        """