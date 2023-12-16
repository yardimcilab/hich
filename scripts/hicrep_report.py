import argparse
import json
import os, subprocess
from pathlib import Path
import matplotlib.pyplot as plt
from ordered_comparisons import *


class ExpSamComparison(Comparison):
    def __init__(self, file1, file2, experiment1, experiment2, filetype):
        if file1 < file2:
            self.files = (file1, file2)
            self.experiments = (experiment1, experiment2)
        else:
            self.files = (file2, file1)
            self.experiments = (experiment2, experiment1)
        self.filetype = filetype
    
    def SameFile(self):
        return self.files[0] == self.files[1]

    def SameExperiment(self):
        return self.experiments[0] == self.experiments[1]

    def __eq__(self, other):
        return self.files == other.files

    def __lt__(self, other):
        if self.SameFile() and not other.SameFile():
            return True
        if other.SameFile() and not self.SameFile():
            return False
        
        if self.filetype == "exp" and other.filetype == "sam":
            return True
        if other.filetype == "exp" and self.filetype == "sam":
            return False
        
        if self.SameExperiment() and not other.SameExperiment():
            return True
        if other.SameExperiment() and not self.SameExperiment():
            return False
        
        if ''.join(self.files) < ''.join(other.files):
            return True
        if ''.join(other.files) < ''.join(self.files):
            return False
        

        return False

    def __str__(self):
        return f"{Path(self.files[0]).stem},\n{Path(self.files[1]).stem}"

    def __hash__(self):
        return hash(self.files)

def save_boxplots(data, binSize, h, dBPMax, filename):
    # Set the figure size
    plt.figure(figsize=(len(data) * .8, 6))  # Adjust the size as needed

    # Prepare the data for plotting
    labels, values = zip(*data)

    # Create the boxplot
    plt.boxplot(values, labels=labels, vert=True)

    # Set the title
    plt.title(f"HiCRep Per-Chromosome SCC Scores\nBin size {binSize}, h {h}, dBPMax {dBPMax}")
    
    # Set x and y axis labels
    plt.xlabel("Comparison")
    plt.ylabel("Per-Chromosome SCC Scores")

    # Set custom x-tick positions and labels
    plt.xticks(ticks=range(1, len(labels) + 1), labels=labels, rotation=45, ha='right')

    # Adjust layout to make room for the x-axis labels and title
    plt.tight_layout()

    # Adjust the bottom margin to ensure x-axis labels are contained in the image
    plt.subplots_adjust(bottom=0.4)

    # Make sure directory exists for file
    Path(filename).parent.mkdir(parents=True, exist_ok=True)

    # Save the plot to a file
    plt.savefig(filename, format='png')

def generate_ordered_comparisons(experiments, compare_to_self):
    comparisons = set()

    exp_pairs = ListPairs(list(experiments.keys()))
    sam_pairs = DictPairs(experiments, experiments)

    for pair in exp_pairs:
        comparisons.add(ExpSamComparison(pair[0], pair[1], pair[0], pair[1], "exp"))
    for pair in sam_pairs:
        comparisons.add(ExpSamComparison(pair[1], pair[3], pair[0], pair[2], "sam"))
    comparisons = sorted(list(comparisons))

    return comparisons

def get_hicrep_scores(filename):
    return [float(line.strip()) for line in open(filename).readlines() if line[0] != "#"]

def generate_hicrep_scores(comparisons, binSizes, h, dBPMax, output_dir):
    score_filenames = []
    data = []

    for binSize in binSizes:
        for comparison in comparisons:
            input_file1, input_file2 = comparison.files
            input_file1 = Path(input_file1)
            input_file2 = Path(input_file2)

            assert input_file1.exists(), f"{input_file1} does not exist"
            assert input_file2.exists(), f"{input_file2} does not exist"
            
            output_filename = Path(f"hicrep/{binSize}/HiCRep_binSize_{binSize}_h_{h}_dBPMax_{dBPMax}_{input_file1.stem}_{input_file2.stem}_.txt")
            if not output_filename.exists():
                output_filename.parent.mkdir(parents=True, exist_ok=True)
                cmd = f"hicrep {input_file1} {input_file2} {output_filename} --binSize {binSize} --h {h} --dBPMax {dBPMax}"
                subprocess.run(cmd, shell=True)

            score_filenames.append(str(output_filename))
            scores = get_hicrep_scores(output_filename)
            data.append((str(comparison), scores))
        save_boxplots(data, binSize, h, dBPMax, Path(output_dir)/str(binSize)/"figures"/f"binSize_{binSize}.png")
        data = []

def main():
    parser = argparse.ArgumentParser(description="Process command line arguments.")

    # Experiments argument: a dict in string format
    parser.add_argument("--experiments", required=True, type=json.loads,
                        help="A dictionary of experiment names with a list of samples.")

    # Experiments argument: a dict in string format
    parser.add_argument("--output_dir", required=True, type=str,
                        help="The directory to store HiCRep results. Figures will be stored in figures output_dir/figures.")

    # binSizes argument: optional, comma-delimited list
    parser.add_argument("--binSizes", type=lambda s: [int(item) for item in s.split(',')],
                        help="Comma-delimited list of bin sizes.", default=[])

    # h argument: integer
    parser.add_argument("--h", required=True, type=int, help="An integer value for 'h'.")

    # dBPMax argument: integer
    parser.add_argument("--dBPMax", required=True, type=int, help="An integer value for 'dBPMax'.")

    # compareToSelf argument: boolean, defaults to True
    parser.add_argument("--compareToSelf", dest="compare_to_self", action="store_false",
                        help="Include comparisons to self. Defaults to True. Set this flag to disable.")

    args = parser.parse_args()

    # Generate comparisons
    comparisons = generate_ordered_comparisons(args.experiments, args.compare_to_self)

    # Generate HiCRep scores
    scores = generate_hicrep_scores(comparisons, args.binSizes, args.h, args.dBPMax, args.output_dir)



if __name__ == "__main__":
    main()
