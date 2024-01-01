include: "common.smk"

configfile: "workflow/config/config.yaml"
MERGES = MergePlan(kv_reverse(config['merge_structure']))


wildcard_constraints:
    read = '|'.join(["1", "2"]),
    replicate = '|'.join(MERGES.replicates()),
    merge = '|'.join(MERGES.merges()),
    downsample = '|'.join([str(ds) for ds in config['downsample']]),
    node1 = '|'.join(MERGES.nodes()),
    node2 = '|'.join(MERGES.nodes())