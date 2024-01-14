rm -rf benchmarks
snakemake --software-deployment-method conda --conda-frontend mamba --cores 24 --until align
