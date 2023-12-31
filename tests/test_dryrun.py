import unittest
import subprocess
from contextlib import contextmanager
import os
from pathlib import Path

@contextmanager
def temporary_cd(path):
    old_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)

class TestDryRun(unittest.TestCase):
    def setUp(self):
        # Set the base directory relative to this script's location
        base_dir = Path(__file__).parent.parent

        # Change directory to base_dir, run Snakemake, then revert directory
        with temporary_cd(base_dir):
            command = ["snakemake", "-n", "-F"]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                raise RuntimeError("Snakemake dry run failed: " + result.stderr)

            # Store the results for use in other test methods
            self.result = result.stdout

    def test_merge(self):
        pass
        #self.assertIn("some expected text", self.result)

if __name__ == '__main__':
    unittest.main()
