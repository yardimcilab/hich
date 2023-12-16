import cooltools, cooler
import hicstraw
import numpy as np
import os
import pandas as pd

bin_sizes = "5000,10000,50000,100000,200000,500000,1000000,3000000,6000000".split(",")

exp = "Mock"
num = "3"
mcool_fn = f"EXP1M_{exp}/merge/matrix/EXP1M_{exp}_ds1.mcool"
hic_fn = mcool_fn.replace(".mcool", ".hic")

data_type = 'observed' # (previous default / "main" data) or 'oe' (observed/expected)
normalization = ["VC", "VC_SQRT", "KR", "SCALE"]  # , VC, VC_SQRT, KR, SCALE, etc.

hic = hicstraw.HiCFile(hic_fn)

for bin_size in bin_sizes:
    assert bin_size in hic.getResolutions(), \
        f"{bin_size} is not part of the possible resolutions {','.join(hic.getResolutions())}"

chrom_sizes = pd.Series({chrom.name: chrom.length for chrom in hic.getChromosomes() if chrom.name != "All"})

# First write the chromosome sizes:
with open(hic.getGenomeID() + '.size', 'w') as fsize:
    for chrom in hic.getChromosomes():
        if chrom.name != "All":
            fsize.write(f"{chrom.name}\t{chrom.length}\n")
# Then write the counts in text file:
with open(cool_file.replace('.cool', ".txt"), 'w') as fo:
    for i in range(len(chrom_sizes)):
        for j in range(i, len(chrom_sizes)):
            chrom1 = chrom_sizes.index[i]
            chrom2 = chrom_sizes.index[j]
            result = hicstraw.straw(data_type, normalization, hic_file, chrom1, chrom2, 'BP', resolution)
            for k in range(len(result)):
                start1 = result[k].binX
                start2 = result[k].binY
                value = result[k].counts
                fo.write(f"{chrom1}\t{start1}\t{start1}\t{chrom2}\t{start2}\t{start2}\t{value}\n")

os.system(f"cooler load -f bg2 {hic.getGenomeID()}.size:{resolution} {cool_file.replace('.cool', '.txt')} {cool_file}")
