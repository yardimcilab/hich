12/17/2023
Ben Skubi, skubi@ohsu.edu
Labs of Gurkan Yardimci and Andrew Adey, OHSU

Snakemake workflow for processing from raw .fastq to per-replicate and merged per-experiment .hic and .mcool with several normalization options (VC, VC_SQRT, KR, SCALE).

Notes on specific rules:

rules/setup.smk:
rule  download_index_reference
    Note: A no-alts reference is essential. See Heng Li's explanation at https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    Note: The reference must be indexed (producing .amb, .ann, .bwt, .pac, .sa). This can take about 3 hours.
    Note: You can either drop symlinks to your reference files in ./resources, or else set parameters to let setup.smk download and index

rules/hic.smk
rule pairtools_parse:
    Note: filter out fragments with mapq < 10 (also in rule hic)
rule pairtools_downsample:
    Note: downsample current replicate .pairs file according to [all replicates min total mapped]/[current replicate total mapped]
rule deduplicate
    Note: remove per-replicate duplicates prior to merging
    Note: this is appropriate for biological replicates; technical replicates should be deduplicated after merging
rule hic 
    Note: generate normalizations and resolutions specified in workflow/config/config.yaml
    Note: throw away reads that map to the same restriction fragment
    Note: filter out fragments with mapq < 10 (also in rule pairtools_parse)
rule mcool: convert .hic to .mcool
    Note (12/17/2023): relies on hic2cool which hasn't been updated to accommodate the latest version of .hic file format. In .hic files, chrom names are just the chromosome number/letter (i.e. 1, 2, X, Y, M, etc). In .mcool files, chrom names are specified per-resolution and should be prefixed with 'chr.' Working with the latest .hic file format, hic2cool prefixes the chromosome numbers at a particular resolution with a variable number of 'chr' repetitions (0-3 have been observed), (i.e. chrchr1, chrchr2). These can be fixed by loading the Cooler, setting the chromnames manually, and saving, but this functionality hasn't been implemented yet. It's unknown if other bugs also result from the conversion.

Time and space requirements:
Processing ~120 Gb of .fastq reads may require about a week with no downsampling.
With 24 cores, it took about 30 minutes to run with two experiments of 3 pairs of 1 million-read .fastq files each

In my testing:
An 868 Mb set of .fastq files generates 47 Gb of intermediate and output files (54x increase).
A 116 Gc set of .fastq files generates 1.8 Tb of intermediate and output files (15x increase).
I haven't tested memory usage.