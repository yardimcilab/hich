import sys
from pathlib import Path

conda:
    "../envs/global.yaml"

def setup_environment():
    # Define the path you want to add to sys.path
    workflow_path = Path(workflow.basedir).parent.as_posix()

    # Add the new path to sys.path if it's not already there
    if workflow_path not in sys.path:
        sys.path.insert(0, workflow_path)


setup_environment()

from scripts.merge_plan import *