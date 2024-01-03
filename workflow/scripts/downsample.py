
import subprocess
from merge_plan import *

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

def compute_replicate_total_mapped(statsfile, replicates):
    total_mapped = {}

    for replicate in replicates:
        filename = statsfile.format(replicate = replicate)
        total_mapped[replicate] = int(read_pairtools_stats(filename, "total_mapped"))

    return total_mapped

def min_downsample(total_mapped, wildcards):
    all_total_mapped = total_mapped.values()
    min_total_mapped = min(all_total_mapped)

    replicate = wildcards.replicate
    replicate_total_mapped = total_mapped[replicate]
    
    return min_total_mapped / replicate_total_mapped

downsample = snakemake.wildcards.downsample
if downsample == "min":
    total_mapped = compute_replicate_total_mapped("results/{replicate}/pairs/{replicate}_pairtools_parse_stats.txt", \
                                                  snakemake.params.replicates.split(" "))
    downsample = min_downsample(total_mapped, snakemake.wildcards)

subprocess.run(f"""pairtools sample \
                    --output \"{snakemake.output}\" \
                    --seed 0 \
                    {downsample} \
                   \"{snakemake.input.parse}\" \
                    &> {snakemake.log}""", shell=True)
