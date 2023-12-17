12/17/2023
Ben Skubi, skubi@ohsu.edu
Labs of Gurkan Yardimci and Andrew Adey, OHSU

Snakemake workflow for processing from raw .fastq to per-replicate and merged per-experiment .hic and .mcool with several normalization options (VC, VC_SQRT, KR, SCALE).

workflow/config/config.yaml: config file for specifying experiment structure and processing parameters
    experiment_structure: associate experiments with biological replicates to control merging
    assembly: hg38 or other reference genome name
    assembly_urls: permits downloading and indexing no alts reference genome automatically 
    chromsizes: path to no alts chromsizes file
    juicer_tools_jar: path to juicer tools jar
    min_mapq: minimum bwa mem mapq filter
    downsample: fraction of .pairs entries to downsample to. 0-1 or min; min results in downsampling to the size of the smallest .pairs replicate
    resolutions: resolutions/bin sizes to produce in the .hic and .mcool outputs

workflow/rules: individual modules containing 1+ snakemake rules; common.smk instead contains common Python functions used by other rules
    setup.smk: downloads resource files
    hic: processes from raw .fastq to .hic and .mcool file formats
    common.smk: contains Python functions and classes used in snakemake rules

resources: reference genomes, chromsizes files, etc

results: output location for intermediate and final files
    final outputs per replicate are in results/{experiment}/{replicate}/matrix



12/14/2023
Ben Skubi, skubi@ohsu.edu
Labs of Gurkan Yardimci and Andrew Adey, OHSU

Snakemake workflow for processing from raw .fastq to per-replicate and merged per-experiment .hic and .mcool with several normalization options (VC, VC_SQRT, KR, SCALE).

********Note 0: Adaptation to SLURM********
The instructions below are how I run this pipeline on an Adey lab node NOT using SLURM.
You'll need to edit this file to be SLURM-compatible if running on Exacloud.

Input:
* paired-end, gzipped fastq files named with a unique per-replicate [prefix], stored/symlinked in [base]/[prefix]_[1 or 2].fq.gz

Processing steps:
* Alignment to reference genome
* Downsampling, MAPQ filtering, deduplication, and selecting 'UU'
* Merging biological replicates

Output:
* Per-replicate and merged .hic and .mcool files with several normalization options
* Intermediate processing files

Each experiment will be contained in a directory named after that experiment.
Per-replicate files will be stored in a subdirectory named after that replicate (.hic and .mcool files will be in [experiment]/[replicate]/matrix/[replicate]_ds[downreplicate].[hic or mcool])
Per-experiment merged files will be stored in a subdirectory named "merge" (.hic and .mcool files will be in [experiment]/[merge]/matrix/[experiment]_ds[downreplicate].[hic or mcool])

Installation:
1. Install dependencies (see below)

To use:
1. Make a base directory [base] containing [base]/hickory.smk, [base]/experiment.py and a subdirectory [base]/fastq/
2. Symlink your fastq files to [base]/fastq; each should be formatted [prefix]_[1 or 2].fq.gz. The original fastqs must be in .gz format.
3. Set EXERIMENTS, REFERENCES, PARAMETERS, and SCRIPTS dicts (see below for details)
4. Run:
    
    snakemake --snakefile hickory.smk --cores 24 --printshellcmds

Feel free to adjust the number of cores

Here is the way hickory.smk is structured:
* EXPERIMENTS is a dict where keys are replicate prefixes and values are associated experiment prefixes. Replicates having the same experiment prefix will be merged.
* REFERENCES contains paths to your reference .fa and reference chromosome .sizes file.
    ******Note that it is ESSENTIAL that you use a "no alts" version of your reference*******
        The reason is that if you have alts, bwa mem will find an enormous number of alternative (XA) alignments to the alt contigs,
        and all these reads will receive a MAPQ of 0 and be filtered out.
* PARAMETERS contains key experimental parameters.
    "assembly" corresponds to the reference you're using.
    "resolutions" is a list of bin sizes you would like to generate for your .hic and .mcool files
    "downreplicate" produces one or more downreplicated versions of your .pairs files, all of which will otherwise be processed identically.
        Use "1" if you want no downsampling.
        Use a number from 0-1 (i.e. "0.5") if you want to downreplicate all .pairs files to the same fraction of their original reads
        Using "min" as the parameter detects MIN_TOTAL_MAPPED, the minimum number of total reads across all per-replicate .pairs files.
            It then computes REPLICATE_i_TOTAL_MAPPED for each replicate REPLICATE_i, then downreplicates REPLICATE_i to MIN_TOTAL_MAPPED/REPLICATE_i_TOTAL_MAPPED.
            To put it simply it more or less downreplicates all per-replicate .pairs to the smallest replicate .pairs file.
* SCRIPTS contains paths to scripts necessary for running hickory.smk. Currently it only requires the juicer_tools .jar file.
    If needed, you can download juicer tools from https://github.com/aidenlab/juicer/wiki/Download

You'll need an environment with the following applications (and their dependencies) installed:
bwa mem
samtools
pairtools
juicer tools
hic2cool

***Note 1: Bug in hic2cool***
hic2cool doesn't handle the newest version of the .hic format perfectly. I don't know the extent of the problem
but it at least misformats the chromosome names. .hic has just the number/name for chromosomes (i.e. 1, 2, 3, ..., X, Y, M).
By contrast, the .mcool produced by hic2cool is supposed to have chromosome names prefixed by 'chr' (i.e. chr1, chr2, ..., chrX).
There are separate chromosome name entries in the .mcool format for each bin size. What I'm seeing is that hic2cool is producing
chromosome name vectors with a variable number of 'chr' prefixes (i.e. 1, 2, 3, or chr1, chr2, chr3, or chrchr1, chrchr2, chrchr3),
with a distinct number of 'chr' prefixes for different bin sizes. As far as I know you can just manually rename them if you're interacting
with a .mcool.

***Note 2: Biological replicates***
This workflow assumes per-replicate replicates are biological replicates; if they are technical replicates, then a per-experiment deduplication step
should ideally be included after merging.

***Note 3: Development***
This workflow is in early development and may contain bugs and other shortcomings. Please contact Ben Skubi with any issues you encounter.

***Note 4: Time and hard drive space requirements***
Processing ~120 Gb of .fastq reads may require about a week with no downsampling.
With 24 cores, it took about 30 minutes to run with two experiments of 3 pairs of 1 million-read .fastq files each

In my testing:
An 868 Mb set of .fastq files generates 47 Gb of intermediate and output files (54x increase).
A 116 Gc set of .fastq files generates 1.8 Tb of intermediate and output files (15x increase).

I haven't tested memory usage.

***Note 5: QC reports***
A second snakefile chicory.smk produces QC reports based on the outputs of hickory.smk, contact Ben for details.
