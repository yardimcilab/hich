import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def load_and_average_scores(filename):
    with open(filename, 'r') as file:
        scores = []
        for line in file:
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    score = float(line)
                    if not np.isnan(score):
                        scores.append(score)
                except ValueError:
                    pass  # Ignore lines that cannot be converted to float
    return np.mean(scores) if scores else np.nan


def parse_caption_file(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 4 and parts not in data:
                data.append(parts)
    return pd.DataFrame(data, columns=['scores_file', 'downsample', 'name1', 'name2'])

def create_heatmap(data, downsample, output_file):
    # Set display options
    pd.set_option('display.max_rows', None)  # or a large number instead of None
    pd.set_option('display.max_columns', None)  # or a large number
    pd.set_option('display.width', None)  # or a specific width
    pd.set_option('display.max_colwidth', None)  # or a specific column width

    duplicates = data.duplicated(subset=['name1', 'name2'])
    if duplicates.any():
        print("Duplicate entries found:", data[duplicates])

    pivot_table = data.pivot(index='name1', columns='name2', values='average_score')
    plt.figure(figsize=(12, 10))  # Adjust figure size as needed
    ax = sns.heatmap(pivot_table, annot=True, fmt=".3g", cmap="cividis", vmin=-1, vmax=1)
    plt.title(f"Mean of Per-Chromosome HiCRep Scores\nDownsample: {downsample}")

    # Rotate column labels and align them
    plt.xticks(rotation=45, ha='right')  # Rotate by 45 degrees and align to the right

    # Ensure the bottom of the heatmap is visible
    plt.tight_layout()

    plt.savefig(output_file)

def main(scores_files, caption_file, output_file):
    scores = {file: load_and_average_scores(file) for file in scores_files}
    captions = parse_caption_file(caption_file)

    captions['average_score'] = captions['scores_file'].map(scores)
    downsample = captions['downsample'].iloc[0]  # Assuming the same downsample for all entries

    print("Num captions:", len(captions))

    create_heatmap(captions, downsample, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate heatmap from HiCRep scores.")
    parser.add_argument('--scores', nargs='+', help='List of score filenames', required=True)
    parser.add_argument('--caption', help='Caption filename', required=True)
    parser.add_argument('--output', help='Output PNG filename', required=True)
    args = parser.parse_args()

    main(args.scores, args.caption, args.output)

