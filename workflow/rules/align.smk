rule align_4DN_hic:
    input:
        r1 = "fastq/{replicate}_1.fq.gz",
        r2 = "fastq/{replicate}_2.fq.gz"
    output:
        "results/{experiment}/sambam/{replicate}.bam"
    shell:
        "bwa mem -SP5M -t {24} {genome} {input.r1} {input.r2} | samtools view -b -o {output}"