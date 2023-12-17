"""Retrieve and format filenames for groups of experimental replicates.

Classes:

    Experiment

Functions:

    generate_kwarg_combinations
    all_subsequences
"""

import itertools as it
from pathlib import Path, PosixPath
import yaml

def kv_reverse(dictionary):
    """Given a dictionary with iterable values, set value elements as keys and associated key as value"""
    reversed = {}
    for key, value in dictionary.items():
        for val in value:
            reversed[val] = key
    return reversed

def generate_kwarg_combinations(multi_kwargs, singleton_kwarg_rule = lambda x: x):
    """Given a dict with 1+ possible kwarg values, return a list of dicts containing all combinations of those values.

    Arguments:
        multi_kwargs (expects dict):
            Keys are keywords
            Values are keyword value options, which can be singletons or collections of singletons.

        singleton_kwarg_rule (default lambda x: x)
            A function that transforms all values according to a user-specifiable rule.
            Useful to define which input values should be treated as singletons vs. collections of singletons.
    
    Returns:
        List of dicts. Each dict has a unique combination of value singletons for each keyword.
    
    Examples:
        multi_kwargs = {"foo":[1, 2] "bar":(3, 4)}
        singleton_kwarg_rule = lambda x: x
        returns: [{"foo":1, "bar":3}, {"foo":1, "bar":4}, {"foo":2, "bar":3}, {"foo":2, "bar":4}]

        multi_kwargs = {"foo":[1, 2] "bar":(3, 4)}
        singleton_kwarg_rule = lambda x: x if type(x) == list else [x]
        returns: [{"foo":1, "bar":(3, 4)}, {"foo":2, "bar":(3, 4)}]
    """
    if multi_kwargs == {}:
        return []
    
    # Extract keys and list of values
    keys, values = zip(*multi_kwargs.items())

    values = [singleton_kwarg_rule(value) for value in values]

    # Generate all combinations of values
    combinations = it.product(*values)

    # Create a list of dictionaries for each combination
    result = [dict(zip(keys, combo)) for combo in combinations]

    return result

def all_subsequences(string):
    """Return all subsequences of the given string"""
    ps = set()
    for s in range(len(string)):
        for e in range(s+1, len(string)+1):
            ps.add(string[s:e])
    return ps

class Experiment(dict):
    """Associate replicates (keys) with experimental conditions (values) and format filenames according to this structure.

    Public methods:
        AssertNameConvention: Assert that all experiment names and their associated replicates have a common, uniquely identifying substring (case-insensitive).
        Replicates: Return a list of replicate names for the given experiments.
        Experiments: Return a list of experiment names for the given replicates.
        FormatTemplate: Generate an iterable for all possible combinations of kwargs and associated experiment/replicate pairs.
    
        Overridden methods:
        items(): Limit or expand which (replicate, experiment) pairs to return based on provided replicate/experiment names and limit/expand keyword.
    
    Public instance variables:
        keys: Replicate names
        values: Experiment name associated with the replicate
    """
    def __init__(self, experiment_structure):
        """Set up the replicate-experiment associtions (experiment_structure contains replicate names as keys, experiment names as values)"""
        super().__init__(experiment_structure)
        self.AssertNameConvention()

    def AssertNameConvention(self):
        """Check that every experiment name and its associated replicate names have a uniquely identifying, case-insensitive substring in common.
        
        Exceptions raised:
            AssertionError if no uniquely identifying substring is found for one or more experiments.
        """
        for experiment in self.Experiments():
            # Generate cis_intersection, the intersection of the substring sets for the names of the experiment and its replicates
            cis_replicates = self.Replicates(experiments = experiment)
            cis = cis_replicates + [experiment]
            cis = [name.upper() for name in cis]
            all_cis_substring_sets = [all_subsequences(name) for name in cis]
            cis_intersection = all_cis_substring_sets[0]
            for cis_substring_set in all_cis_substring_sets:
                cis_intersection = cis_intersection.intersection(cis_substring_set)

            # Generate trans_union, the union of the substring sets for the names of other experiments and their replicates
            other_experiments = self.Experiments().difference(set([experiment]))
            trans_replicates = self.Replicates(experiments = other_experiments)
            trans = trans_replicates + list(other_experiments)
            trans = [name.upper() for name in trans]
            all_trans_substring_sets = [all_subsequences(name) for name in trans]
            trans_union = set()
            for trans_substring_set in all_trans_substring_sets:
                trans_union = trans_union.union(trans_substring_set)
            
            # Substrings uniquely identifying the current experiment and its replicates exist if the difference between
            # cis_intersection and trans_union is not empty.
            unique_substrings = cis_intersection.difference(trans_union)
            assert unique_substrings != set(), "No substring specific to experiment and its replicates found in prefixes for {experiment}! You may have a typo in your prefixes."

    def Replicates(self, experiments = None):
        """Returns a list of replicate names associated with the given experiment names (or all replicate names if experiments is None)"""
        if experiments is None:
            experiments = self.values()
        return [replicate for replicate, experiment in super().items() if experiment in experiments]
    
    def Experiments(self, replicates = None):
        """Returns the set of all experiment names associated with the given replicate names (or all experiment names if replicates is None)"""
        if replicates is None:
            replicates = self.keys()
        return set([self[replicate] for replicate in replicates])
    
    def items(self, replicate = "all", experiment = "all", type = "limit"):
        """Returns list of tuples of (replicate, experiment) name pairs.
        
        Arguments:
            replicate: A single replicate name or list of them (defaults to "all").

            experiment: A single experiment name or list of them (defaults to "all").

            type (default "limit"):
                If "limit" is specified, returns only (replicate, experiment) pairs where both the replicate and its experiment are provided
                If "expand" is specified, returns every (replicate, experiment) pair where either the replicate or its experiment are provided
        """
        result = []
        if replicate == "all":
            replicate = self.Replicates()
        if experiment == "all":
            experiment = self.Experiments()
        if type == "limit":
            result = [(rep, exp) for rep, exp in super().items() if rep in replicate and exp in experiment]
        elif type == "expand":
            result = [(rep, exp) for rep, exp in super().items() if rep in replicate or exp in experiment]
        return result
    
    def FormatTemplate(self, templates, singleton_kwarg_rule = lambda x: x if type(x) in [list, set] else [x], rep_exp_relation = "limit", **kwargs):
        """Format template(s) with stored {experiment}, {replicate} and other kwargs, and return an iterable of (formatted_string, replicate, experiment) tuples.
        
        Arguments:
            template: a stringifiable object or list/set/tuple of stringifiable objects that may contain wildcards.

            singleton_kwarg_rule: enclose kwarg value in list if it is a singleton. (default: anything other that list or set is a singleton)

            rep_exp_relation: if "limit", only expands {replicate} if its experiment is provided. If "expand", expands {replicate} with all possibilities for provided values of {experiment}.

            **kwargs: keys are keywords, values are singletons (converted into single-element collections) or collections of singletons.
        
        Returns:
            Iterable tuple (formatted_string, replicate, experiment). Formats with all replicates and experiments by default. If experiments provided in **kwargs,
            formats with just the provided experiments and their associated samples. If samples provided in **kwargs, will format with just those samples. If both
            are provided in **kwargs, will only format with samples in **kwargs that also have their associated experiments in **kwargs. Generates all possible combinations
            of provided **kwargs except in the case of experiment and replicate names, where replicate names are only provided with their associated experiment. Applies
            singleton_kwarg_rule on kwarg values to decide if values are singletons or collections of singletons for purposes of forming combinations.
        """
        # If a bare string is submitted, put it in a list for iterating in the main loop
        
        if type(templates) not in [list, set]:
            templates = [templates]

        for template in templates:
            # Convert to string (i.e. to handle pathlib Paths and PosixPaths)
            template_string = str(template)
            
            # If replicate or experiment provided in kwargs, separate them out from the rest of the kwargs
            replicate = "all" if kwargs.get('replicate') is None else kwargs.get('replicate')
            experiment = "all" if kwargs.get('experiment') is None else kwargs.get('experiment')
            kwargs.pop('replicate', None)
            kwargs.pop('experiment', None)
            
            # Generate combinations of non-replicate/experiment kwargs and of replicate/experiment combinations to use
            combinations = generate_kwarg_combinations(kwargs, singleton_kwarg_rule)
            items = self.items(replicate = replicate, experiment = experiment, type = rep_exp_relation)

            for replicate, experiment in items:
                base_kwargs = {'replicate': replicate, 'experiment': experiment}

                # There's probably a way to simplify this...
                if combinations:
                    for unique_kwargs in combinations:
                        # Update the non-replicate/experiment unique keyword combination with the current replicate/experiment values
                        unique_kwargs.update(base_kwargs)

                        formatted_string = template_string.format(**unique_kwargs)
                        yield (formatted_string, replicate, experiment)
                else:
                    formatted_string = template_string.format(**base_kwargs)
                    yield (formatted_string, replicate, experiment)

EXPERIMENTS = Experiment(kv_reverse(config['experiment_structure']))


wildcard_constraints:
    replicate = '|'.join(EXPERIMENTS.Replicates()),
    experiment = '|'.join(EXPERIMENTS.Experiments()),
    downsample = '|'.join([str(ds) for ds in config['downsample']])

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

def compute_replicate_total_mapped(statsfile):
    pairtools_parse_stats = EXPERIMENTS.FormatTemplate(statsfile)
    total_mapped = {}

    for filename, replicate, experiment in pairtools_parse_stats:
        total_mapped[replicate] = int(read_pairtools_stats(filename, "total_mapped"))

    return total_mapped

def min_downreplicate(total_mapped, wildcards):
    all_total_mapped = total_mapped.values()
    min_total_mapped = min(all_total_mapped)

    replicate = wildcards.replicate
    replicate_total_mapped = total_mapped[replicate]
    
    return min_total_mapped / replicate_total_mapped

"""
exp = Experiment({"E1S1":"E1", "E1S2":"E1", "E2S1":"E2", "E2S2":"E2"})
exp.AssertNameConvention()
print(list(exp.FormatTemplate("{foo} {bar}", foo = [1, 2], bar = (3, 4))))
print(list(exp.FormatTemplate("{foo} {bar}", foo = [1, 2], bar = (3, 4), singleton_kwarg_rule=lambda x: x )))
"""

"""
print(generate_kwarg_combinations({"foo":[1, 2], "bar":(3, 4)}))
print(generate_kwarg_combinations({"foo":[1, 2], "bar":(3, 4)}, lambda x: x if type(x) == list else [x]))
"""

"""
EXPERIMENTS = Experiment({"1M_KO_1"  :"EXP1M_KO",
                          "1M_KO_2"  :"EXP1M_KO",
                          "1M_KO_3"  :"EXP1M_KO",
                          "1M_Mock_1":"EXP1M_Mock",
                          "1M_Mock_2":"EXP1M_Mock",
                          "1M_Mock_3":"EXP1M_Mock"})

EXPERIMENTS.AssertNameConvention()
print(EXPERIMENTS.ReplicateExperimentPairs(replicate = "1M_KO_1", experiment = "EXP1M_Mock", type = "limit"))
print(EXPERIMENTS.ReplicateExperimentPairs(replicate = "1M_KO_1", experiment = "EXP1M_Mock", type = "expand"))
result = list(EXPERIMENTS.FormatTemplate("Experiment: {experiment} Replicate: {replicate} a: {a} b: {b}", experiment="EXP1M_Mock", replicate=["1M_KO_1", "1M_KO_2"], a=[1, 2], b=[3, 4], rep_exp_relation="expand"))

for r in result:
    print(r)
"""