"""Retrieve and format filenames for groups of experimental replicates.

Classes:

    Merge

Functions:

    generate_kwarg_combinations
    all_subsequences
"""

from networkx import DiGraph
import itertools as it
from pathlib import Path, PosixPath
import yaml

def kv_reverse(dictionary):
    """Given a dictionary with iterable values, set value elements as keys and associated key as value"""
    reversed = {}
    for key, values in dictionary.items():
        for val in values:
            reversed.setdefault(val, set())
            reversed[val].add(key)
    return reversed

def kwarg_product(**kwargs):
    """Given kwargs with values as iterables,
    return a list of dicts with the same keys
    and all combinations of individual values"""

    """
    Create lists of (key, value) tuples for each key and each individual value
    Get the power set of tuples
    Set them as key, value parts in a dict
    Save a list of these dicts
    Return that list of dicts
    """
    kwargs_options = []
    for key, values in kwargs.items():
        if len(values) > 0:
            options = [(key, value) for value in values]
            kwargs_options.append(options)
    return set(it.product(*kwargs_options))

class MergePlan(DiGraph):
    def __init__(self, merge_structure):
        super().__init__()
        """Set up the replicate-merge associtions (merge_structure contains replicate names as keys, merge names as values)"""
        for key, values in merge_structure.items():
            self.add_node(key)
            for value in values:
                self.add_node(value)
                self.add_edge(key, value)
    
    def format_template_group_combinations(self, template, groups, combo_kwarg_keys, with_replacement = True, crash_toosmall_groups = True, **kwargs):
        formatted = []
        comb_func = it.combinations_with_replacement if with_replacement else it.combinations
        r = len(combo_kwarg_keys)
        for group in groups:
            if len(group) < r:
                assert not crash_toosmall_groups, f"{len(group)} elements in {group}, but combination size is set at {r} by number of kwarg_keys {combo_kwarg_keys}!"
                continue
            combinations = sorted(list(comb_func(group, r)))
            for combination in combinations:
                new_kwargs = {key: [val] for key, val in zip(combo_kwarg_keys, sorted(combination))}
                new_kwargs.update(kwargs)
                formatted += list(self.format_template(template, **new_kwargs))
        return sorted(formatted)

    def format_template(self, template, **kwargs):
        rep_mrg_options = self._rep_mrg_options(**kwargs)
        result = set()
        for rmo in rep_mrg_options:
            rep_mrg = rmo[0]
            options = rmo[1]
            fmt = dict(options)
            fmt["replicate"] = rep_mrg[0]
            fmt["merge"] = rep_mrg[1]
            result.add(template.format(**fmt))
        return result
        
    def _rep_mrg_options(self, **kwargs):
        rep_mrg, other_kwargs = self._rep_mrg(**kwargs)
        other_kwarg_product = kwarg_product(**other_kwargs)
        rep_mrg_options = it.product(rep_mrg.edges(), other_kwarg_product)
        return rep_mrg_options
    
    def _rep_mrg(self, **kwargs):
        kwargs_replicates, kwargs_merges = self._get_rep_mrg_or_empty(**kwargs)

        replicates = self.replicates(kwargs_merges)
        merges = self.merges(kwargs_replicates)

        if kwargs_replicates != set():
            replicates = replicates.intersection(kwargs_replicates)
        if kwargs_merges != set():
            merges = merges.intersection(kwargs_merges)

        kwargs.pop("replicates", None)
        kwargs.pop("merge", None)
        
        return self._rep_mrg_intersection(replicates, merges), kwargs

    def _get_rep_mrg_or_empty(self, **kwargs):
        kwarg_replicates = set(kwargs.get("replicates", set()))
        kwarg_merges = set(kwargs.get("merge", set()))
        return kwarg_replicates, kwarg_merges

    def _rep_mrg_intersection(self, replicates, merges):
        G = {}
        for replicate in replicates:
            keep_merges = set(self.successors(replicate)).intersection(merges)
            if len(keep_merges) > 0:
                G[replicate] = keep_merges
        return MergePlan(G)

    def replicates(self, merges = set()):
        reps = set()
        if merges == set():
            return self.source_nodes()
        for merge in merges:
            reps = reps.union(set(self.predecessors(merge)))
        return reps
    
    def merges(self, replicates = set()):
        mrgs = set()
        if replicates == set():
            return self.sink_nodes()
        for replicate in replicates:
            mrgs = mrgs.union(set(self.successors(replicate)))
        return mrgs
    
    def same_structure(self, other):
        return set(self.nodes()) == set(other.nodes()) and set(self.edges()) == set(other.edges())
    
    def source_nodes(self):
        return set([node for node in self.nodes() if self.out_degree(node) > 0 and self.in_degree(node) == 0])
    
    def sink_nodes(self):
        return set([node for node in self.nodes() if self.in_degree(node) > 0 and self.out_degree(node) == 0])