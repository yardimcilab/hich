from HiCkory.workflow.scripts.rep_merge import *
import unittest

class TestMergeStructure(unittest.TestCase):
    def setUp(self):
        self.merge = Merge({
            "a1": "A",
            "a2": "A",
            "a3": "A",
            "b1": "B",
            "b2": "B",
            "b3": "B"
        })

    def test_all_replicates(self):
        assert set(self.merge.replicates()) == set(["a1", "a2", "a3", "b1", "b2", "b3"])
    
    def test_subset_replicates(self):
        assert set(self.merge.replicates(["B"])) == set(["b1", "b2", "b3"])
    
    def test_all_merges(self):
        assert set(self.merge.merges()) == set(["A", "B"])
    
    def test_subset_replicates_uniform(self):
        assert set(self.merge.merges(["a1", "a2"])) == set(["A"])
    
    def test_subset_replicates_both(self):
        assert set(self.merge.merges(["a1", "b2"])) == set(["A", "B"])
    
    def test_format_template(self):
        fmt = set(self.merge.format_template("{merge}{replicate}"))
        assert fmt == set([("Aa1", "a1", "A"), \
                            ("Aa2", "a2", "A"), \
                            ("Aa3", "a3", "A"), \
                            ("Bb1", "b1", "B"), \
                            ("Bb2", "b2", "B"), \
                            ("Bb3", "b3", "B")])