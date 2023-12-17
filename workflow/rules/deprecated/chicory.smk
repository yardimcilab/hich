from pathlib import Path
import itertools as it
from statistics import mean
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cooltools, cooler

import warnings
warnings.filterwarnings("ignore")

# import semi-core packages
from matplotlib import colors
import numpy as np
import pandas as pd
plt.style.use('seaborn-poster')

import bioframe

EXPERIMENTS = {
    "EXP1M_KO": ["1M_KO_1", "1M_KO_2", "1M_KO_3"],
    "EXP1M_Mock": ["1M_Mock_1", "1M_Mock_2", "1M_Mock_3"]
}

PROCESSING = {
    "downsample": [1],
    "bin_size": "5000,10000,50000,100000,200000,500000,1000000,3000000,6000000".split(","),
    "hicrep_h": 1,
    "hicrep_dbpmax": 6000000
}

REFERENCES = {
    "genome": "/home/groups/ravnica/refs/hg38/hg38.fa",
    "chromsizes": "/home/groups/ravnica/projects/sciMET_GCC/230410_replication/ref/hg38.noalt.sizes"
}

DIR = {
    "fastq": Path("fastq"),
    "sample_matrix": Path("{experiment}") / "{sample}" / "matrix",
    "experiment_matrix": Path("{experiment}") / "merge" / "matrix",
    "fastqc": Path("chicory") / "fastqc",
    "decay_curve_figures": Path("chicory") / "decay_curve" / "ds{downsample}" / "bin_size{bin_size}",
    "hicrep_scores": Path("chicory") / "hicrep" / "ds{downsample}" / "bin_size{bin_size}",
    "hicrep_figures": Path("chicory") / "hicrep" / "ds{downsample}" / "bin_size{bin_size}" / "figures"
}

INPUT_FILES = {
    "fastq":                    [DIR["fastq"] / "{sample}_1.fq.gz", DIR["fastq"] /  "{sample}_2.fq.gz"],
    "mcool":                    [DIR["sample_matrix"] / "{sample}_ds{downsample}.mcool"],
    "mcool_merge":              [DIR["experiment_matrix"] / "{experiment}_ds{downsample}.mcool"]
}

OUTPUT_FILES = {
    "fastqc":                   [DIR["fastqc"] / "{sample}_1_fastqc.html", DIR["fastqc"] / "{sample}_2_fastqc.html"],
    "decay_curve_figures":      [DIR["decay_curve_figures"] / "DecayCurve_{prefix}.png"],
    "hicrep_scores":            [DIR["hicrep_scores"] / "HiCRep_{prefix1}_{prefix2}_ds{downsample}_bin_size{bin_size}.txt"],
    "hicrep_figures":           [DIR["hicrep_figures"] / "HiCRep_ds{downsample}_bin_size{bin_size}.png"]
}

SCRIPTS = {
    "fastqc": "/home/groups/ravnica/src/fastqc/FastQC-0.11.9/fastqc"
}

# Use bioframe to fetch the genomic features from the UCSC.
hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_chromsizes = pd.DataFrame({
    'chrom': hg38_chromsizes.index,
    'length': hg38_chromsizes.values
})
hg38_cens = bioframe.fetch_centromeres('hg38')
# create a view with chromosome arms using chromosome sizes and definition of centromeres
hg38_arms = bioframe.make_chromarms(hg38_chromsizes,  hg38_cens)

# Generate combinations of experiments

def find_experiment(sample):
    for exp, sams in EXPERIMENTS.items():
        if sample in sams:
            return exp
    assert False, (f"Could not find experiment associated with sample {sample}")

ALL_PATHS = []
exp_exp = list(it.combinations_with_replacement(EXPERIMENTS.keys(), 2))
sam_sam = list(it.combinations_with_replacement(it.chain.from_iterable(EXPERIMENTS.values()), 2))
comparisons = exp_exp + sam_sam
for downsample in PROCESSING["downsample"]:
    for exp, samps in EXPERIMENTS.items():
        ALL_PATHS += expand(list(it.chain.from_iterable(INPUT_FILES.values())), experiment=[exp], sample=samps, downsample=downsample)
        ALL_PATHS += expand(OUTPUT_FILES["fastqc"], sample=samps)
        ALL_PATHS += expand(OUTPUT_FILES["decay_curve_figures"], prefix=[exp]+samps, downsample=PROCESSING['downsample'], bin_size=PROCESSING['bin_size'])
    for hicrep_comparison in comparisons:
        ALL_PATHS += expand(OUTPUT_FILES["hicrep_scores"], prefix1=hicrep_comparison[0], prefix2=hicrep_comparison[1], downsample=downsample, bin_size=PROCESSING["bin_size"])
        ALL_PATHS += expand(OUTPUT_FILES["hicrep_figures"], downsample=downsample, bin_size=PROCESSING["bin_size"])
        

wildcard_constraints:
    sample = '|'.join(it.chain.from_iterable(EXPERIMENTS.values())),
    experiment = '|'.join(EXPERIMENTS.keys()),
    prefix = '|'.join(list(it.chain.from_iterable(EXPERIMENTS.values())) + list(EXPERIMENTS.keys())),
    prefix1 = '|'.join(list(it.chain.from_iterable(EXPERIMENTS.values())) + list(EXPERIMENTS.keys())),
    prefix2 = '|'.join(list(it.chain.from_iterable(EXPERIMENTS.values())) + list(EXPERIMENTS.keys())),
    bin_size = '|'.join(PROCESSING["bin_size"])
    

rule all:
    input:
        ALL_PATHS

rule fastqc:
    input:
        INPUT_FILES["fastq"]
    output:
        OUTPUT_FILES["fastqc"]
    run:
        output_path = Path(output[0]).parent
        output_path.mkdir(parents=True, exist_ok=True)
        for fastq in input:
            shell(f"{SCRIPTS['fastqc']} -o {output_path} {fastq}")

def get_decay_curve_figure_inputs(wildcards):
    if wildcards.prefix in EXPERIMENTS.keys():
        return str(INPUT_FILES["mcool_merge"][0]).format(experiment=wildcards.prefix, downsample=wildcards.downsample)
    else:
        return str(INPUT_FILES["mcool"][0]).format(experiment=find_experiment(wildcards.prefix), sample=wildcards.prefix, downsample=wildcards.downsample)

def get_decay_curve_figure_outputs(wildcards):
    return str(OUTPUT_FILES["decay_curve_figures"][0]).format(prefix=wildcards.prefix, downsample=wildcards.downsample, bin_size=PROCESSING["bin_size"])

rule decay_curve_figures:
    input:
        lambda wildcards: get_decay_curve_figure_inputs(wildcards)
    output:
        [str(OUTPUT_FILES["decay_curve_figures"][0]).format(bin_size=bin_size, downsample="{downsample}", prefix="{prefix}") for bin_size in PROCESSING["bin_size"]]
    run:
        # Code adapted from cooltools tutorial on decay curves
        # https://cooltools.readthedocs.io/en/latest/notebooks/contacts_vs_distance.html
        
        for bin_size, output_filename in zip(PROCESSING["bin_size"], output):
            clr = cooler.Cooler(f"{input}::/resolutions/{bin_size}")
            Path(output[0]).parent.mkdir(exist_ok=True, parents=True)

            # hic2cool appears to screw up chromosome names, so we fix here
            rename_chrs = {}
            
            for chrID in clr.chromnames:
                new_chrID = "chr" + chrID.replace("chr", "")
                rename_chrs[chrID] = new_chrID

            cooler.rename_chroms(clr, rename_chrs)

            # select only those chromosomes available in cooler
            arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

            # cvd == contacts-vs-distance
            cvd = cooltools.expected_cis(
                clr=clr,
                view_df=arms,
                smooth=False,
                aggregate_smoothed=False,
                clr_weight_name='SCALE',
                nproc=1 #if you do not have multiple cores available, set to 1
            )

            cvd['s_bp'] = cvd['dist']* int(bin_size)
            f, ax = plt.subplots(1,1)

            for region in hg38_arms['name']:
                ax.loglog(
                    cvd['s_bp'].loc[cvd['region1']==region],
                    cvd['balanced.avg'].loc[cvd['region1']==region],
                )
                ax.set(
                    xlabel='separation, bp',
                    ylabel='IC contact frequency')
                ax.set_aspect(1.0)
                ax.grid(lw=0.5)
            plt.savefig(str(output_filename))

def hicrep_scores_get_input(wildcards):
    if wildcards.prefix1 in EXPERIMENTS.keys():
        s1 = str(INPUT_FILES["mcool_merge"][0]).format(experiment=wildcards.prefix1, downsample=wildcards.downsample, bin_size=wildcards.bin_size)
        s2 = str(INPUT_FILES["mcool_merge"][0]).format(experiment=wildcards.prefix2, downsample=wildcards.downsample, bin_size=wildcards.bin_size)
        return (s1, s2)
    elif wildcards.prefix1 in it.chain.from_iterable(EXPERIMENTS.values()):
        s1 = str(INPUT_FILES["mcool"][0]).format(experiment = find_experiment(wildcards.prefix1), sample=wildcards.prefix1, downsample=wildcards.downsample, bin_size=wildcards.bin_size)
        s2 = str(INPUT_FILES["mcool"][0]).format(experiment = find_experiment(wildcards.prefix2), sample=wildcards.prefix2, downsample=wildcards.downsample, bin_size=wildcards.bin_size)
        return (s1, s2)
    assert False, f"Problem with wildcards: {wildcards}"

rule hicrep_scores:
    input:
        lambda wildcards: hicrep_scores_get_input(wildcards)
    output:
        OUTPUT_FILES["hicrep_scores"]
    run:
        cmd = f"hicrep {input} {output} --binSize {wildcards.bin_size} --h {PROCESSING['hicrep_h']} --dBPMax {PROCESSING['hicrep_dbpmax']}"
        shell(cmd)

def hicrep_figures_get_input(wildcards):
    return [str(OUTPUT_FILES["hicrep_scores"][0]).format(prefix1=c[0], prefix2=c[1], downsample=wildcards.downsample, bin_size=wildcards.bin_size) for c in comparisons]

rule hicrep_figures:
    input:
        lambda wildcards: hicrep_figures_get_input(wildcards)
    output:
        OUTPUT_FILES["hicrep_figures"]
    run:
        data = {}
        for comparison, scc_file in zip(comparisons, input):
            scc_score = mean([float(line.strip()) for line in open(scc_file).readlines() if line[0] != '#'])
            data[comparison] = scc_score
        # Extract unique elements for y and x axes
        y_labels = sorted(set(key[0] for key in data.keys()))
        x_labels = sorted(set(key[1] for key in data.keys()))

        # Create an empty DataFrame
        df = pd.DataFrame(index=y_labels, columns=x_labels)

        # Fill the DataFrame with your data
        for (y, x), value in data.items():
            df.at[y, x] = value

        # Create the heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(df.astype(float), annot=True, cmap='coolwarm')
        plt.title(f"Heatmap of HiCRep Per-Chromosome SCC Score Means\nBin size: {wildcards.bin_size} H: {PROCESSING['hicrep_h']} dBPMax: {PROCESSING['hicrep_dbpmax']} downsample: {wildcards.downsample}")
        plt.ylabel("Matrix 1")
        plt.xlabel("Matrix 2")

        # Adjust the bottom margin
        plt.subplots_adjust(bottom=0.2)  # Increase the bottom margin size


        # Save the heatmap as a PNG file
        plt.savefig(str(output))
        plt.close()