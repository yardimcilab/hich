rm -rf benchmarks
snakemake --software-deployment-method conda --conda-frontend mamba "$@"
