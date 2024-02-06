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

# Rule all

HiCkory's outputs are governed by the files you tell it to generate in `rule all` in your Snakefile.

Producing the per-condition `.hic` and `.mcool` files will also force generation of all preprocessing files, but not the per-technical replicate or per-biological replicate `.hic` or `.mcool` files. These must be explicitly required as inputs to `rule all`.

You must also explicitly require reports to be generated, such as `fastqc` or `multiqc`.

In the future, a more convenient interface will be created for selecting which outputs should be generated.

# Workflow Structure and Config Specification

Adapting the example `config.yaml` file for your experiment should be straightforward. Here is a detailed description of what HiCkory expects its structure to be.

## data:
The data section currently just has one key, `fastq:`, which specifies the relative path where HiCkory will look for input .fastq files.

**Formatting fastq file names**

Pay attention to the `conditions:`, `bioreps:`, and `techreps:` sections below. HiCkory expects fastq files to be given names corresponding to the names specified in these sections according to the format:

`[base HiCkory directory]/[fastq]/[cond]_[biorep]_[techrep]_1.fq.gz`

`[base HiCkory directory]/[fastq]/[cond]_[biorep]_[techrep]_2.fq.gz`

**Example**

If you have a condition named `cond1` with a biological replicate named `br1` with a technical replicate named `tr1`, and have specified `fastq: "fastq/"`, then you must name the paired-end fastq files `fastq/cond1_br1_tr1_1.fq.gz` and `fastq/cond1_br1_tr1_2.fq.gz`. Note that `.gz` indicates the fastq files have been compressed with [gzip](https://www.gzip.org/).

See `workflow/rules/techreps.smk` in `rule bwamem2_bam` to see exactly how HiCkory looks for the input fastq files based on these settings.

## processes:
Each process represents a single set of HiCkory parameters.
If you gather more data and wish to add new replicates or conditions to the workflow, or you wish to rerun with different parameters, you can define these changes in `config.yaml` within the `processes:` list.

Each process has two required sections: `parameters` and `conditions`. You can also specify other sections *ad hoc* to specify metadata or whatever else you like.

### parameters:
The parameters section describes processing parameters that are applied across all conditions.

#### min_mapq
Minimum MAPping Quality (MAPQ) thresholds applied at different rules. In principle, MAPQ describes the likelihood of a read being incorrectly aligned as `MAPQ = -10log10(ChanceOfWrongAlignment)`, such that mapq = 60 represents a 1/1000000 chance of an incorrect alignment. In practice, each aligner provides an approximation of MAPQ using an aligner-specific range. BWA MEM(2) has a range from 0-60.

##### align
MAPQ threshold to output alignments from the aligner.
Note: in future versions that permit a choice of aligner, this will be moved to a per-condition parameter.

##### parse
MAPQ threshold to retain pairs during `rule parse`.

#### downsample
How to downsample reads. There are two choices:
##### to_constant_fraction:
Downsample to a constant fraction. If `1` or unspecified, size will be unchanged.
###### for_techreps: `[0-1]`
###### for_bioreps: `[0-1]`
###### for_conditions: `[0-1]`

##### to_min_pairs_size:
After defining a comparison group, downsample to the size of the smallest sample in the comparison group. Comparison groups can be specified for technical replicates, biological replicates, and conditions.

The calculation that uses these parameters is performed in `workflow/scripts/calculate_downsample.py`.

For the given level of processing, .pairs files will be compared only within the grouping specified. For example, if performing `to_min_pairs_size` downsampling on biological replicates using `compare: within_conditions`, each biological replicate files will be downsampled to the size of the smallest biological replicate .pairs file within the same condition.

###### for_techreps:
`  - compare: [within_processes|within_conditions|within_bioreps|within_techreps]`
###### for_bioreps:
`  - compare: [within_processes|within_conditions|within_bioreps|within_techreps]`
###### for_conditions:
`  - compare: [within_processes|within_conditions|within_bioreps|within_techreps]`


#### matrix
Parameters governing creation of .hic and .mcool files
##### bin_sizes
List of resolutions/bin sizes to create. Should be a list of integers specified as strings (i.e. in quotation marks).
##### digest
True/False depending on whether you wish to produce an additional 'resolution' that is per-restriction fragment rather than a uniform bin size. If True, you must provide an appropriately formatted restriction fragment file in the per-condition `genome/restriction_sites` section.
##### norms
Which contact matrix norms to compute in the `.hic` contact matrices in addition to raw counts.
### conditions:
The conditions section specifies the structure of each experimental condition, including condition-specific processing parameters and identifiers for technical and biological replicates.

In addition to the required parameters described here, you can add additional *ad hoc* parameters to store information such as notes.
#### `[condition name]`
The name of the condition. Should be unique.
##### genome:
Information on how to align the condition.
###### index:
Path to the index used for alignment. This typically should be a no-alts index. See Heng Li's commentary at https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use

Note that BWA MEM and BWA MEM2 use different aligners. HiCkory uses BWA MEM2, so you'll need to create a BWA MEM2 index.
###### assembly:
Keyword for genome assembly used. Parameter for `pairtools parse` where it can be anything and is just a piece of metadata. It is a functional parameter in `juicer_tools pre`, see https://github.com/aidenlab/juicer/wiki/Pre for a list of options. If your option isn't listed, you can use the path to the chrom.sizes file.
###### chromsizes:
Path to chromosome sizes file
###### restriction_sites:
`  - juicer_tools:` Path to restriction sites file formatted for juicer tools
Some kits (i.e. Arima) provide a script for generating this file.

Format:

`chromosome1 cutsite1 cutsite2 cutsite 3...`

`chromosome2 cutsite1 cutsite2 cutsite 3...`

Example:

`chr10 10444 10550 10693 10867 10962 12128...`

`chr10_GL383545v1_alt 44 49 428 1114 1202 1397...`

`  - cooltools:` Path to restriction sites file formatted for cooltools

Format:
`chromosome start end name 0 +`

Example:
`chr1	0	11159	HIC_chr1_1	0	+`

##### bioreps:
List of biological replicate names
###### techreps:
List of technical replicate names