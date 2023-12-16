from pathlib import Path
import json, copy, subprocess
import itertools as it
from experiment import Experiment

#***********************************
# Customization section
#***********************************

# Replace list with unique prefixes for your biological replicates and keys with unique prefixes for merging per-experiment replicates
# All fastq files should be in gzip format with the filename [prefix]_[1 or 2].fq.gz and stored/symlinked at [base]/fastq/[prefix]_[1 or 2].fq.gz
# By default, each experiment prefix needs a substring unique to it and its corresponding replicates (case insensitive; KO and Mock here)
# This is to help ensure you don't merge the wrong replicates by accident by associating replicates and experiments incorrectly.
EXPERIMENTS = Experiment({"1M_KO_1"  :"EXP1M_KO",
                          "1M_KO_2"  :"EXP1M_KO",
                          "1M_KO_3"  :"EXP1M_KO",
                          "1M_Mock_1":"EXP1M_Mock",
                          "1M_Mock_2":"EXP1M_Mock",
                          "1M_Mock_3":"EXP1M_Mock"})

# Replace with filenames for reference genome (NO ALTS!!!) and chromsizes file
REFERENCES = {
    "genome": "/home/groups/ravnica/refs/hg38/hg38.fa",
    "chromsizes": "/home/groups/ravnica/projects/sciMET_GCC/230410_replication/ref/hg38.noalt.sizes"
}

# Replace with assembly name for your reference, the size of the contact matrix resolutions you want to produce, downsampling strategies, and minimum mapq to filter for
# downreplicate "1" = no downsampling, "min" = downreplicate each replicate by (min total mapped reads across all replicates)/(total mapped reads in current replicate)
# can also set a value of "0.5" or any other fraction to downreplicate by a uniform fraction
PARAMETERS = {
    "assembly": "hg38",
    "resolutions": [1E3, 5E3, 1E4, 5E4, 1E5, 2E5, 5E5, 1E6, 3E6, 6E6],
    "downreplicate": [1, "min"],
    "min_mapq": 10
}

# Replace with path to your juicer_tools_[version].jar file
# Download and installation link: https://github.com/aidenlab/juicer/wiki/Download
SCRIPTS = {
    "juicer_tools_pre": "java -Xmx20g -jar /home/groups/ravnica/src/juicer/juicer_tools_1.22.01.jar pre"
}

#***********************************
# Static section
#***********************************
# Do not alter below unless you know what you're doing, this defines the file structure for the inputs and outputs
# Nothing below this point should need customization.
DIR = {
    "fastq": Path("fastq"),
    "reports": Path("{experiment}") / "{replicate}" / "reports",
    "sambam": Path("{experiment}") / "{replicate}" / "sambam",
    "replicate_pairs": Path("{experiment}") / "{replicate}" / "pairs",
    "replicate_matrix": Path("{experiment}") / "{replicate}" / "matrix",
    "experiment_pairs": Path("{experiment}") / "merge" / "pairs",
    "experiment_matrix": Path("{experiment}") / "merge" / "matrix"
}

FILES = {
    "fastq":                    [DIR["fastq"] / "{replicate}_1.fq.gz", DIR["fastq"] /  "{replicate}_2.fq.gz"],
    "align":                    [DIR["sambam"] / "{replicate}.bam"],
    "name_sort":                [DIR["sambam"] / "{replicate}.name_sort.bam"],
    "pairtools_parse":          [DIR["replicate_pairs"] / "{replicate}.pairs"],
    "pairtools_parse_stats":    [DIR["replicate_pairs"] / "{replicate}_pairtools_parse_stats.txt"],
    "pairtools_replicate":         [DIR["replicate_pairs"] / "{replicate}_ds{downreplicate}.pairs"],
    "pairtools_sort":           [DIR["replicate_pairs"] / "{replicate}_ds{downreplicate}_sort.pairs"],
    "pairtools_deduplicate":    [DIR["replicate_pairs"] / "{replicate}_ds{downreplicate}_sort_deduplicate.pairs"],
    "pairtools_select":         [DIR["replicate_pairs"] / "{replicate}_ds{downreplicate}_sort_deduplicate_select.pairs"],
    "hic":                      [DIR["replicate_matrix"] / "{replicate}_ds{downreplicate}.hic"],
    "mcool":                    [DIR["replicate_matrix"] / "{replicate}_ds{downreplicate}.mcool"],
    "pairtools_merge":          [DIR["experiment_pairs"] / "{experiment}_ds{downreplicate}_merge.pairs"],
    "pairtools_sort_merge":     [DIR["experiment_pairs"] / "{experiment}_ds{downreplicate}_merge_sort.pairs"],
    "pairtools_select_merge":   [DIR["experiment_pairs"] / "{experiment}_ds{downreplicate}_merge_sort_select.pairs"],
    "hic_merge":                [DIR["experiment_matrix"] / "{experiment}_ds{downreplicate}.hic"],
    "mcool_merge":              [DIR["experiment_matrix"] / "{experiment}_ds{downreplicate}.mcool"]
}

EXPERIMENTS.AssertNameConvention()

wildcard_constraints:
    replicate = '|'.join(EXPERIMENTS.Replicates()),
    experiment = '|'.join(EXPERIMENTS.Experiments()),
    downreplicate = '|'.join([str(ds) for ds in PARAMETERS["downreplicate"]])

def AllPaths():
    generic_paths = it.chain.from_iterable(FILES.values())
    unique_generic_paths = set(generic_paths)
    
    downreplicate = PARAMETERS["downreplicate"]

    formatted_paths = EXPERIMENTS.FormatGeneric(unique_generic_paths, downreplicate = downreplicate)
    
    all_paths = [filename for filename, _, _ in formatted_paths]
    return all_paths

rule all:
    input:
        AllPaths()

rule align:
    input:
        FILES["fastq"]
    output:
        FILES["align"]
    log:
        "log/align/align{experiment}_{replicate}.log"
    run:
        genome = REFERENCES["genome"]
        r1, r2 = input
        threads = 24
        params_4DN = "-SP5M"
        convert_to_bam = "samtools view -b"

        shell(f"bwa mem {params_4DN} -t {threads} {genome} {r1} {r2} | {convert_to_bam} -o {output}")


rule name_sort:
    input:
        rules.align.output
    output:
        FILES["name_sort"]
    run:
        shell(f"samtools sort -n -o --deliberateerror eroijog {output} {input}")

def read_pairtools_stats(filename, query_name):
    all_lines = open(filename).readlines()

    for line in all_lines:
        split = line.split()

        if len(split) < 2:
            continue
        
        line_name = split[0].strip()
        line_stat = split[1].strip()

        if query_name == line_name:
            return line_stat

rule pairtools_parse:
    input:
        rules.name_sort.output
    output:
        pairtools_parse = FILES["pairtools_parse"],
        pairtools_parse_stats = FILES["pairtools_parse_stats"]
    run:
        chromsizes = REFERENCES['chromsizes']
        assembly = PARAMETERS['assembly']
        shell(f"""
            pairtools parse {input} \
                -o {output.pairtools_parse} \
                --output-stats {output.pairtools_parse_stats} \
                -c {chromsizes} \
                --assembly {assembly} \
                --min-mapq {PARAMETERS["min_mapq"]} \
                --nproc-in 24 \
                --nproc-out 24
            """)

# To refactor - this is ugly
rule wait_for_output:
    input:
        [filename for filename, _, _ in EXPERIMENTS.FormatGeneric(FILES["pairtools_parse"][0])]
    output:
        marker = "all_replicates_pairtools_parsed"
    run:
        shell(f"touch {output.marker}")

def compute_replicate_total_mapped():
    pairtools_parse_stats = EXPERIMENTS.FormatGeneric(FILES["pairtools_parse_stats"][0])
    PARAMETERS.setdefault("replicate_total_mapped", {})

    for filename, replicate, experiment in pairtools_parse_stats:
        total_mapped = int(read_pairtools_stats(filename, "total_mapped"))
        PARAMETERS["replicate_total_mapped"][replicate] = total_mapped

def min_downreplicate(wildcards):
    all_total_mapped = PARAMETERS["replicate_total_mapped"].values()
    min_total_mapped = min(all_total_mapped)

    replicate = wildcards.replicate
    replicate_total_mapped = PARAMETERS["replicate_total_mapped"][replicate]
    
    return min_total_mapped / replicate_total_mapped

rule pairtools_replicate:
    input:
        files = rules.pairtools_parse.output.pairtools_parse,
        marker = rules.wait_for_output.output
    output:
        FILES["pairtools_replicate"]
    run:
        compute_replicate_total_mapped()

        pairtools_parse_generic = FILES["pairtools_parse"][0]
        input_filename, _, _ = next(EXPERIMENTS.FormatGeneric(pairtools_parse_generic, **wildcards))

        downreplicate = wildcards.downreplicate if wildcards.downreplicate != "min" else min_downreplicate(wildcards)
        
        shell(f"pairtools replicate --output {output} --seed 0 {downreplicate} {input_filename}")

rule pairtools_sort:
    input:
        FILES["pairtools_replicate"]
    output:
        FILES["pairtools_sort"]
    run:
        shell(f"pairtools sort -o {output} {input}")

rule pairtools_deduplicate:
    input:
        rules.pairtools_sort.output
    output:
        FILES["pairtools_deduplicate"]
    run:
        shell(f"pairtools dedup --mark-dups -o {output} {input}")

def format_selection(pair_types):
    specifiers = []

    for pair_type in pair_types:
        specifier = f"(pair_type==\"{pair_type}\")"
        specifiers.append(specifier)

    select = ' or '.join(specifiers)

    return f"'{select}'"

rule pairtools_select:
    input:
        rules.pairtools_deduplicate.output
    output:
        FILES["pairtools_select"]
    run:
        select = format_selection(["UU", "UR", "RU"])
        shell(f"pairtools select {select} -o {output} {input}")

def int2str(num_list, sep=','):
    strlist = [str(int(value)) for value in num_list]
    return sep.join(strlist)

rule hic:
    input:
        rules.pairtools_select.output
    output:
        FILES["hic"]
    run:
        juicer_tools_pre = SCRIPTS["juicer_tools_pre"]
        resolutions = int2str(PARAMETERS["resolutions"])
        assembly = PARAMETERS['assembly']

        shell(f"{juicer_tools_pre} -r {resolutions} {input} {output} {assembly}")

rule mcool:
    input:
        rules.hic.output
    output:
        FILES["mcool"]
    run:
        shell(f"hic2cool convert {input} {output} -p 24")

def pairtools_select_to_merge(wildcards):
    generic_filename = FILES["pairtools_select"]
    experiment = wildcards.experiment
    downreplicate = wildcards.downreplicate
    
    formatted = EXPERIMENTS.FormatGeneric(generic_filename, experiment = experiment, downreplicate = downreplicate)

    return [filename for filename, _, _ in formatted]

rule pairtools_merge:
    input:
        lambda wildcards: pairtools_select_to_merge(wildcards)
    output:
        FILES["pairtools_merge"]
    run:
        for i in input:
            print(i)
        shell(f"pairtools merge -o {output} {input}")

rule pairtools_sort_merge:
    input:
        rules.pairtools_merge.output
    output:
        FILES["pairtools_sort_merge"]
    run:
        shell(f"pairtools sort -o {output} {input}")

rule select_merge:
    input:
        rules.pairtools_sort_merge.output
    output:
        FILES["pairtools_select_merge"]
    run:
        select = format_selection(["UU", "UR", "RU"])
        shell(f"pairtools select {select} -o {output} {input}")

rule hic_merge:
    input:
        rules.select_merge.output
    output:
        FILES["hic_merge"]
    run:
        juicer_tools_pre = SCRIPTS["juicer_tools_pre"]
        resolutions = int2str(PARAMETERS["resolutions"])
        assembly = PARAMETERS['assembly']
        
        merge_file_dir = Path(output[0]).parent
        merge_file_dir.mkdir(exist_ok=True, parents=True)

        shell(f"{juicer_tools_pre} -r {resolutions} {input} {output} {assembly}")

rule mcool_merge:
    input:
        rules.hic_merge.output
    output:
        FILES["mcool_merge"]
    run:
        shell(f"hic2cool convert {input} {output} -p 24")

    