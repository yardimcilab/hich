from HiCkory.workflow.scripts.rep_merge import *
import unittest

class TestMergePlan(unittest.TestCase):
    def setUp(self):
        self.merge = MergePlan({
            "a1": ["A"],
            "a2": ["A", "A23"],
            "a3": ["A", "A23"],
            "b1": ["B"],
            "b2": ["B"],
            "b3": ["B"]
        })

        self.simple_merge = MergePlan({
            "a1": ["A"],
            "a2": ["A"],
            "b1": ["B"]
        })
    
    def test_format_template(self):
        assert self.simple_merge.format_template("{merge}{replicate}{opt1}", opt1 = ["o1", "o2"]) == \
            set(["Aa1o1", "Aa1o2", "Aa2o1", "Aa2o2", "Bb1o1", "Bb1o2"])
        assert self.simple_merge.format_template("{merge}{replicate}", opt1 = ["o1", "o2"]) == \
            set(["Aa1", "Aa2", "Bb1"])

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
        rep_mrg, kwargs = self.merge._rep_mrg(replicates = ["a1"], merges = ["A"])
        assert rep_mrg.same_structure(MergePlan({"a1": ["A"]}))
        assert kwargs == {}

        rep_mrg, kwargs = self.merge._rep_mrg(replicates = ["a1", "a2"], merges = ["A", "A23"], moil = ["for gold"])
        assert rep_mrg.same_structure(MergePlan({"a1": ["A"], "a2": ["A", "A23"]}))
        assert kwargs == {"moil": ["for gold"]}

        rep_mrg, kwargs = self.merge._rep_mrg()
        assert rep_mrg.same_structure(self.merge)
        assert kwargs == {}

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