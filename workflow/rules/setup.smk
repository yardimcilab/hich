rule download_juicer_tools:
    output:
        "resources/juicer_tools"
    shell:
        """wget -O resources/juicer_tools https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar"""
