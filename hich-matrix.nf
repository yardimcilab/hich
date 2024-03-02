include { CAT } from './hickory-common.nf'

params.juicer_tools = file('resources/juicer_tools')
params.bins = [1000000]
params.assembly = "hg38"
params.chromsizes = file('resources/hg38_noalts.sizes')
params.results = 'results'

/*
https://github.com/aidenlab/juicer/wiki/Pre
*/
process HIC_MATRIX {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'

    input:
    tuple \
        val(sample_id), \
        path(pairs)
    
    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.hic")
    
    script:
    """
    lz4 -c ${pairs} > ${sample_id}.pairs
    java -Xmx20g -jar \
    ${params.juicer_tools} pre \
        -r ${params.bins.join(",")} \
        ${sample_id}.pairs \
        ${sample_id}.hic \
        ${params.assembly}
    """
}

process COOL_MATRIX {
    tag "${sample_id}"
    publishDir params.results, mode:'copy'

    input:
    tuple \
        val(sample_id), \
        path(pairs)
    
    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.cool")
    
    script:
    """
    lz4 ${pairs} > ${sample_id}.pairs
    cooler cload pairs \
        --assembly ${params.assembly} \
        -c1 2 -p1 3 -c2 4 -p2 5 \
        ${params.chromsizes}:${params.bins.min()} \
        ${sample_id}.pairs \
        ${"${sample_id}.cool"}
    """
}

process MCOOL_MATRIX {
    tag "${sample_id}"
    publishDir params.results, model:'copy'

    input:
    tuple \
        val(sample_id), \
        path(cool)
    
    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.mcool")
    
    script:
    """
    cooler zoomify -r ${params.bins.join(",")} --balance --balance-args '--max-iters 1000' ${cool}
    """
}

/*
The creation of a single-cell cooler file is similar to a regular cooler file. Each cell needs to have a name, bin table and a pixel table. All cells must have the same dimensions, and the bins and pixels needs to be provided as two dicts with the cell names as keys.

sparse 2D arrays, identical axes in HDF5

Bin table: chrom start end [, weight]
    order: chrom (enum), start

Elements/pixels table: bin1_id bin2_id value
Used by HiC-Pro

Also includes chromosome description table and indexes to support random access
chroms: name length

weights are for normalization/matrix balancing, can have multiple or bin-level masks

You can access cooler.Cooler.chroms(), cooler.Cooler.bins(), cooler.Cooler.pixels() from cool files
So one way to make a .scool file should be to make Coolers from the pairs in python,
then access their bins and pixels.

name_pixel_dict = {'cell1': pixels_cell1, 'cell2': pixels_cell2, 'cell3': pixels_cell3}
name_bins_dict = {'cell1': bins_cell1, 'cell2': bins_cell2, 'cell3': bins_cell3}
cooler.create_scool('single_cell_cool.scool', name_bins_dict, name_pixel_dict)
*/

process SCOOL_MATRIX {
    publishDir params.results, model:'copy'

    input:
        tuple \
            val(scool_id), \
            val(block_id), \
            path(cools)
    
    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    import cooler

    scool_file = "${scool_id}.scool"
    block_ids = "${block_id.join(" ")}".split(" ")
    cools = "${cools}".split(" ")
    
    bins = {}
    pixels = {}
    for id, cool in zip(block_ids, cools):
        c = cooler.Cooler(f"{cool}::/")
        bins[id] = c.bins()
        pixels[id] = c.pixels()
    cooler.create_scool(scool_file, bins, pixels)

    """
}

workflow {
    csv_ch = CAT(Channel.fromPath('matrix.csv')).splitCsv().map {tuple(it[0], file(it[1]))}
    
    
    HIC_MATRIX(csv_ch)
    cool_ch = COOL_MATRIX(csv_ch)
    MCOOL_MATRIX(cool_ch)

    scool_input_ch = Channel.from("all").combine(cool_ch).groupTuple()

    scool_ch = SCOOL_MATRIX(scool_input_ch)
    scool_ch.view()
}   