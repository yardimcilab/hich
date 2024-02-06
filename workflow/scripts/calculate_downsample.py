import subprocess

wc = snakemake.wildcards

def writelog(message):
    subprocess.run(f"""echo "{message}" > {snakemake.log}""", shell=True)

def downsample(output, pair_filenames, params):
    subprocess.run(f"""echo "{output} {pair_filenames} {params}" > {snakemake.log}""", shell=True)
    process = params.process
    group_for = f'for_{params.input_type}'
    if 'downsample' not in process['parameters']:
        downsample_to_constant_fraction(output, pair_filenames, 1)
    else:
        ds = process['parameters']['downsample']
        assert not ('to_constant_fraction' in ds and 'to_min_pairs_size' in ds), \
        "Must choose to_constant_fraction or to_min_pairs_size in downsampling parameter, not both"
        
        if 'to_constant_fraction' in ds:
            if group_for not in ds['to_constant_fraction']:
                downsample_to_constant_fraction(output, pair_filenames, 1)
            else:
                downsample_to_constant_fraction(output, pair_filenames, ds['to_constant_fraction'][group_for])
        elif 'to_min_pairs_size' in ds:
            if group_for not in ds['to_min_pairs_size']:
                downsample_to_constant_fraction(output, pair_filenames, 1)
            else:
                downsample_to_min_pairs_size(output, pair_filenames, process, group_for)

def downsample_to_constant_fraction(output, pair_filenames, fraction):
    contents = [f"{file}\t{float(fraction)}" for file in pair_filenames]
    contents_str = '\n'.join(contents)
    open(str(output), "w").write(contents_str)

def downsample_to_min_pairs_size(output, pair_filenames, process, group_for):
    compare = process['parameters']['downsample']['to_min_pairs_size'][group_for]['compare']
    writelog(f"Comparing {compare}")
    filename_groups = group_by_comparison_type(pair_filenames, compare)
    writelog(f"Filename groups {filename_groups}")
    pair_sizes = compute_pair_sizes(pair_filenames)
    writelog(f"Pair sizes {pair_sizes}")
    downsample_fractions = compute_downsample_fractions(filename_groups, pair_sizes)
    writelog(f"Downsample fractions {downsample_fractions}")
    write_downsample_fractions(output, downsample_fractions)

def group_by_comparison_type(pair_filenames, compare):
    compare_to_depth = {
        'within_processes':2,
        'within_conditions':3,
        'within_bioreps':4,
        'within_techreps':5
    }
    depth = compare_to_depth[compare]
    min_depth = min([f.count('/') for f in pair_filenames])
    assert min_depth >= depth, f'During downsampling, minimum file path depth was only {min_depth}, but {depth} is required to compare by {compare}.'
    comparisons = {}
    for f in pair_filenames:
        s = f.split('/')
        start = '/'.join(s[:depth])
        comparisons.setdefault(start, [])
        comparisons[start].append(f)
    return comparisons

def pair_size(filename):
    start = 0
    end = 0
    with open(filename) as f:
        line_number = 0
        for line in f:
            line_number += 1
            if line.startswith("#columns"):
                start = line_number
        if start != 0:
            end = line_number - start
        else:
            raise ValueError(f"#columns not found in {filename}")
    return end

def compute_pair_sizes(pair_filenames):
    pair_sizes = {}
    for f in pair_filenames:
        writelog(f"Computing pair sizes for {f}")
        pair_sizes[f] = pair_size(f)
    return pair_sizes

def compute_downsample_fractions(filename_groups, pair_sizes):
    downsample_fractions = {}
    for group, files in filename_groups.items():
        min_size = pair_sizes[files[0]]
        for file in files:
            min_size = min(min_size, pair_sizes[file])
        for file in files:
            fraction = min_size/pair_sizes[file]
            downsample_fractions[file] = fraction
    return downsample_fractions

def write_downsample_fractions(output, downsample_fractions):
    with open(str(output), "w") as f:
        data = []
        for filename, fraction in downsample_fractions.items():
            data.append(f"{filename}\t{fraction}")
        f.write('\n'.join(data))

downsample(snakemake.output, snakemake.input, snakemake.params)
