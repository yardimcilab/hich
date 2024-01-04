import sys
import os

# Get the absolute path to the directory containing merge_plan.py
script_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'workflow', 'scripts'))

# Add this directory to sys.path
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Now you can import merge_plan
from merge_plan import *


import unittest

class TestMergePlan(unittest.TestCase):
    def setUp(self):
        self.merge_rev_dict = {"A": ["a1", "a2", "a3"], "A23": ["a2", "a3"], "B": ["b1", "b2", "b3"]}
        self.merge_dict = {
            "a1": {"A"},
            "a2": {"A", "A23"},
            "a3": {"A", "A23"},
            "b1": {"B"},
            "b2": {"B"},
            "b3": {"B"}
        }

        self.merge = MergePlan(self.merge_dict)

        self.simple_merge = MergePlan({
            "a1": ["A"],
            "a2": ["A"],
            "b1": ["B"]
        })
    
    def test_format_template_group_combinations(self):
        assert self.simple_merge.format_template_group_combinations("{node1}{node2}{param}", \
                                                                    [self.simple_merge.source_nodes(), self.simple_merge.sink_nodes()], \
                                                                    ["node1", "node2"], \
                                                                    False, \
                                                                    param = ["param"]) == \
                ["ABparam", "a1a2param", "a1b1param", "a2b1param"]
        assert self.simple_merge.format_template_group_combinations("{node1}{node2}", \
                                                                    [self.simple_merge.source_nodes(), self.simple_merge.sink_nodes()], \
                                                                    ["node1", "node2"], \
                                                                    True) == \
                ["AA", "AB", "BB", "a1a1", "a1a2", "a1b1", "a2a2", "a2b1", "b1b1"]
        

    def test_kv_reverse(self):
        assert kv_reverse(self.merge_rev_dict) == self.merge_dict
    
    def test_format_template(self):
        assert self.simple_merge.format_template("{merge}{replicate}{opt1}", opt1 = ["o1", "o2"]) == \
            set(["Aa1o1", "Aa1o2", "Aa2o1", "Aa2o2", "Bb1o1", "Bb1o2"])
        assert self.simple_merge.format_template("{merge}{replicate}", opt1 = ["o1", "o2"]) == \
            set(["Aa1", "Aa2", "Bb1"])
        assert self.merge.format_template("{merge}{replicate}") == \
            set(["Aa1", "Aa2", "Aa3", "A23a2", "A23a3", "Bb1", "Bb2", "Bb3"])
        #assert self.merge.format_template("{replicate}", merge = ["A23"]) == \
        #    set(["a2", "a3"])

    def test_rep_mrg_options(self):
        assert set(self.simple_merge._rep_mrg_options(opt1 = ["o1", "o2"])) == \
            set([(("a1", "A"), (("opt1", "o1"),)), \
                (("a1", "A"), (("opt1", "o2"),)), \
                (("a2", "A"), (("opt1", "o1"),)), \
                (("a2", "A"), (("opt1", "o2"),)), \
                (("b1", "B"), (("opt1", "o1"),)), \
                (("b1", "B"), (("opt1", "o2"),))])

    def test_kwarg_product(self):
        assert kwarg_product(arg1 = ["a1a", "a1b"], arg2 = ["b1a"]) == \
            set([(("arg1", "a1a"), ("arg2", "b1a")), (("arg1", "a1b"), ("arg2", "b1a"))])
        assert kwarg_product(arg1 = [], arg2 = []) == set([()])
        assert kwarg_product(arg1 = ["a1a", "a1b"], arg2 = []) == \
            set([(("arg1", "a1a"), ), (("arg1", "a1b"),)])
    
    def test_rep_mrg(self):
        rep_mrg, kwargs = self.merge._rep_mrg(replicates = ["a1"], merge = ["A"])
        assert rep_mrg.same_structure(MergePlan({"a1": ["A"]}))
        assert kwargs == {}

        rep_mrg, kwargs = self.merge._rep_mrg(replicates = ["a1", "a2"], merge = ["A", "A23"], moil = ["for gold"])
        assert rep_mrg.same_structure(MergePlan({"a1": ["A"], "a2": ["A", "A23"]}))
        assert kwargs == {"moil": ["for gold"]}

        rep_mrg, kwargs = self.merge._rep_mrg()
        assert rep_mrg.same_structure(self.merge)
        assert kwargs == {}

        rep_mrg, kwargs = self.merge._rep_mrg(merge = ["A23"])
        assert rep_mrg.same_structure(MergePlan({"a2": ["A23"], "a3": ["A23"]}))
        assert kwargs == {}
    
    def test_rep_mrg_or_empty(self):
        assert self.merge._get_rep_mrg_or_empty(**{"merge": ["A23"]}) == (set(), {"A23"})

    def test_rep_mrg_intersection(self):
        assert self.merge._rep_mrg_intersection(["a1"], ["A"]).same_structure(MergePlan({"a1": ["A"]}))
        assert self.merge._rep_mrg_intersection(["a1", "a2"], ["A"]).same_structure( \
            MergePlan({"a1": ["A"], "a2": ["A"]}))
        assert self.merge._rep_mrg_intersection(["a1", "a2"], ["A23"]).same_structure( \
            MergePlan({"a2": ["A23"]}))
        assert self.merge._rep_mrg_intersection(["a1", "a2"], ["A", "A23"]).same_structure( \
            MergePlan({"a1": ["A"], "a2": ["A", "A23"]}))
        assert self.merge._rep_mrg_intersection(["a1", "a2"], ["A", "A23", "B"]).same_structure( \
            MergePlan({"a1": ["A"], "a2": ["A", "A23"]}))
        assert self.merge._rep_mrg_intersection(["a1", "a2"], []).same_structure( \
            MergePlan({}))
        assert self.merge._rep_mrg_intersection([], ["A", "A23", "B"]).same_structure( \
            MergePlan({}))
        assert self.merge._rep_mrg_intersection([], ["in", "the", "midnight", "sun"]).same_structure( \
            MergePlan({}))
        assert self.merge._rep_mrg_intersection(["a1", "a2", "a3", "b1", "b2", "b3"], ["A23"]).same_structure( \
            MergePlan({"a2": ["A23"], "a3": ["A23"]}))

    def test_replicates(self):
        assert self.merge.replicates() == set(["a1", "a2", "a3", "b1", "b2", "b3"])
        assert self.merge.replicates([]) == set([])
        assert self.merge.replicates(["A23"]) == set(["a2", "a3"])
        assert self.merge.replicates(["A", "A23"]) == set(["a1", "a2", "a3"])
        assert self.merge.replicates(["B", "A23"]) == set(["a2", "a3", "b1", "b2", "b3"])
    
    def test_merges(self):
        assert self.merge.merges() == set(["A", "A23", "B"])
        assert self.merge.merges(["a1"]) == set(["A"])
        assert self.merge.merges(["a2"]) == set(["A", "A23"])
        assert self.merge.merges(["a1", "b1"]) == set(["A", "B"])
        assert self.merge.merges(["a2", "b2"]) == set(["A", "A23", "B"])
    
    def test_same_structure(self):
        assert self.merge.same_structure(self.merge)
        assert not self.merge.same_structure(MergePlan({
            "a1": ["A", "A23"],
            "a2": ["A", "A23"],
            "a3": ["A", "A23"],
            "b1": ["B"],
            "b2": ["B"],
            "b3": ["B"]
        }))
        assert not self.merge.same_structure(MergePlan({}))

    def test_source_nodes(self):
        assert self.merge.source_nodes() == set(["a1", "a2", "a3", "b1", "b2", "b3"])
    
    def test_sink_nodes(self):
        assert self.merge.sink_nodes() == set(["A", "A23", "B"])