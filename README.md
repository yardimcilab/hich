HiCkory is a Hi-C data processing pipeline.

Features
- Based in Snakemake
- Fastq -> .mcool and .hic contact matrices
- Merge technical and biological replicates for multiple conditions
- Downsample to minimum replicate sample size
- Makes fastqc and Hi-C specific multiqc reports
- Dependency management built in with mamba and docker

Future features
- HiCRep clustermaps
- Bisulfite-converted DNA alignment and methylation calling for methyl-HiC
- Post-processing analysis (compartment, loop and TAD calling; promoter-enhancer interactions)
- Config generation
- Pseudo-bulking single-cell data
