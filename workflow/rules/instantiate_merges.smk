include: "common.smk"

MERGES = MergePlan(kv_reverse(config['merge_structure']))

wildcard_constraints:
    replicate = '|'.join(MERGES.replicates()),
    merge = '|'.join(MERGES.merges()),
    downsample = '|'.join([str(ds) for ds in config['downsample']])