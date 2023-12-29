"""Retrieve and format filenames for groups of experimental replicates.

Classes:

    Merge

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

def all_subsequences(string):
    """Return all subsequences of the given string"""
    ps = set()
    for s in range(len(string)):
        for e in range(s+1, len(string)+1):
            ps.add(string[s:e])
    return ps

class Merge(dict):
    """Associate replicates (keys) with experimental conditions (values) and format filenames according to this structure.

    Public methods:
        AssertNameConvention: Assert that all merge names and their associated replicates have a common, uniquely identifying substring (case-insensitive).
        Replicates: Return a list of replicate names for the given merges.
        Experiments: Return a list of merge names for the given replicates.
        FormatTemplate: Generate an iterable for all possible combinations of kwargs and associated merge/replicate pairs.
    
        Overridden methods:
        items(): Limit or expand which (replicate, merge) pairs to return based on provided replicate/merge names and limit/expand keyword.
    
    Public instance variables:
        keys: Replicate names
        values: Merge name associated with the replicate
    """
    def __init__(self, experiment_structure):
        """Set up the replicate-merge associtions (experiment_structure contains replicate names as keys, merge names as values)"""
        super().__init__(experiment_structure)

    def replicates(self, merges = None):
        """Returns a list of replicate names associated with the given merge names (or all replicate names if merges is None)"""
        if merges is None:
            merges = self.values()
        return [replicate for replicate, merge in super().items() if merge in merges]
    
    def merges(self, replicates = None):
        """Returns the set of all merge names associated with the given replicate names (or all merge names if replicates is None)"""
        if replicates is None:
            replicates = self.keys()
        return set([self[replicate] for replicate in replicates])
    
    def items(self, replicate = "all", merge = "all", type = "limit"):
        """Returns list of tuples of (replicate, merge) name pairs.
        
        Arguments:
            replicate: A single replicate name or list of them (defaults to "all").

            merge: A single merge name or list of them (defaults to "all").

            type (default "limit"):
                If "limit" is specified, returns only (replicate, merge) pairs where both the replicate and its merge are provided
                If "expand" is specified, returns every (replicate, merge) pair where either the replicate or its merge are provided
        """
        result = []
        if replicate == "all":
            replicate = self.replicates()
        if merge == "all":
            merge = self.merges()
        if type == "limit":
            result = [(rep, exp) for rep, exp in super().items() if rep in replicate and exp in merge]
        elif type == "expand":
            result = [(rep, exp) for rep, exp in super().items() if rep in replicate or exp in merge]
        return result
    

    
    def format_template(self,
                        templates,
                        singleton_kwarg_rule = lambda x: x if type(x) in [list, set] else [x],
                        rep_exp_relation = "limit",
                        **kwargs):
        """Format template(s) with stored {merge}, {replicate} and other kwargs,
        and return an iterable of (formatted_string, replicate, merge) tuples.
        
        Arguments:
            template: a stringifiable object or list/set/tuple
            of stringifiable objects that may contain wildcards.

            singleton_kwarg_rule: enclose kwarg value in list if it is a singleton.
            (default: anything other that list or set is a singleton)

            rep_exp_relation: if "limit", only expands {replicate} if its merge is provided.
            If "expand", expands {replicate} with all possibilities for provided values of {merge}.

            **kwargs: keys are keywords, values are singletons (converted into single-element collections)
            or collections of singletons.
        
        Returns:
            Iterable tuple (formatted_string, replicate, merge).
            Formats with all replicates and merges by default.
            If merges provided in **kwargs, formats with just
            the provided merges and their associated samples.
            If samples provided in **kwargs, will format with just those samples.
            If both are provided in **kwargs, will only format with samples
            in **kwargs that also have their associated merges in **kwargs.
            Generates all possible combinations of provided **kwargs except
            in the case of merge and replicate names, where replicate names
            are only provided with their associated merge. Applies
            singleton_kwarg_rule on kwarg values to decide if values are singletons
            or collections of singletons for purposes of forming combinations.
        """
        templates = self._ensure_iterable(templates)
        replicate, merge, other_kwargs = self._extract_special_kwargs(kwargs)
        combinations = self._generate_kwarg_combinations(other_kwargs, singleton_kwarg_rule)
        items = self.items(replicate=replicate, merge=merge, type=rep_exp_relation)

        return self._format_and_yield(templates, combinations, items)

    def _generate_kwarg_combinations(self, multi_kwargs, singleton_kwarg_rule = lambda x: x):
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

    def _ensure_iterable(self, templates):
        return templates if isinstance(templates, (list, set)) else [templates]

    def _extract_special_kwargs(self, kwargs):
        replicate = kwargs.pop('replicate', 'all')
        merge = kwargs.pop('merge', 'all')
        return replicate, merge, kwargs

    def _format_and_yield(self, templates, combinations, items):
        for template in templates:
            for replicate, merge in items:
                for unique_kwargs in combinations or [{}]:
                    yield self._format_tuple(template, replicate, merge, unique_kwargs)

    def _format_tuple(self, template, replicate, merge, kwargs):
        base_kwargs = {'replicate': replicate, 'merge': merge}
        kwargs.update(base_kwargs)
        formatted_string = str(template).format(**kwargs)
        return formatted_string, replicate, merge
    
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

def compute_replicate_total_mapped(statsfile, merges):
    pairtools_parse_stats = merges.format_template(statsfile)
    total_mapped = {}

    for filename, replicate, merge in pairtools_parse_stats:
        total_mapped[replicate] = int(read_pairtools_stats(filename, "total_mapped"))

    return total_mapped

def min_downsample(total_mapped, wildcards):
    all_total_mapped = total_mapped.values()
    min_total_mapped = min(all_total_mapped)

    replicate = wildcards.replicate
    replicate_total_mapped = total_mapped[replicate]
    
    return min_total_mapped / replicate_total_mapped