import subprocess

def get_fractions(file):
    fractions = {}
    with open(file) as f:
        for line in f.readlines():
            s = line.split()
            fractions[s[0]] = float(s[1])    
    return fractions

def downsample(output, input, fraction, log):
    cmd = f"""pairtools sample \
                        --output \"{output}\" \
                        --seed 0 \
                        {fraction} \
                        \"{input}\" \
                        1> {log} 2> {log};"""
    subprocess.run(cmd, shell=True)

fractions = get_fractions(snakemake.input.fractions)
downsample(snakemake.output, snakemake.input.pairs, fractions[snakemake.input.pairs], snakemake.log)