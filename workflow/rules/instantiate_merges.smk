include: "common.smk"

MERGES = Merge(kv_reverse(config['merge_structure']))

wildcard_constraints:
    replicate = '|'.join(MERGES.replicates()),
    merge = '|'.join(MERGES.merges()),
    downsample = '|'.join([str(ds) for ds in config['downsample']])