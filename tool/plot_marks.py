import argparse
from pathlib import Path
import subprocess
import struct

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    # Parsing args
    parser = argparse.ArgumentParser(
        description="Plot marks density",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    script_dir = Path(__file__).resolve().parent
    parser.add_argument("-c", "--cmd_path", default=script_dir / "build", help="Benchmarks command path")
    parser.add_argument("-o", "--output_path", default="./", help="Output path")
    parser.add_argument("-g", "--group", action="store_true", help="Group of collections")
    parser.add_argument("collection", help="Collection path")
    args = parser.parse_args()

    data = {}
    collection_path = Path(args.collection)
    if args.group:
        # Processing collections group
        for child in collection_path.iterdir():
            if child.is_dir():
                process_collection(child, args.output_path, args.cmd_path, data)
    else:
        # Processing single collection
        process_collection(collection_path, args.output_path, args.cmd_path, data)

    plot_marks_densities(data, args.output_path)


def esc(code):
    return f'\033[{code}m'


def process_collection(collection_path, output_path, bm_cmd_path, data):
    collection_path = Path(collection_path).resolve()
    output_path = Path(output_path).resolve()
    bm_cmd_path = Path(bm_cmd_path).resolve()

    # Creating output directory for the given collection
    collection_name = collection_path.name
    print(f"{esc('38;5;120')}Collection '{collection_name}'{esc(0)}")
    output_path = output_path / collection_name
    output_path.mkdir(parents=True, exist_ok=True)

    data_path = collection_path / "sri" / "bwt_run_first_text_pos_data.sdsl.vec.bin"
    if not data_path.exists():
        print(f"{esc('38;5;69')}Create file for marks vector{esc(0)}")

        cmd = str(bm_cmd_path / "int_vector_to_vector")
        cmd += " --data=./data"
        cmd += " --key=bwt_run_first_text_pos"
        cmd += " 2>int_vector_to_vector-error.txt"

        subprocess.run(cmd, shell=True, cwd=collection_path / "sri")

    values = read_data(data_path)
    # values = np.array(read_data(data_path)).astype(np.float64)
    values.sort()
    # values *= 100.0 / values[len(values) - 1]

    data[collection_name] = values

    # plot_marks_density(values, collection_name, output_path)


def read_data(data_path):
    values = []
    with open(data_path, 'rb') as fin:
        for chunk in iter(lambda: fin.read(8), b''):
            values.append(struct.unpack('Q', chunk)[0])

    return values


def plot_marks_density(values, collection_name, output_path):
    plot_path = output_path / f"{collection_name}-marks-density"

    for bw_adjust in [0.25, 0.5, 0.75, 1]:
        sns.kdeplot(values, fill=True, bw_adjust=bw_adjust)
        plt.savefig(f"{plot_path}-kde-{bw_adjust}.png")
        plt.clf()

    sns.distplot(a=values, bins=int(len(values) / 100), hist=True, kde=True, rug=False)
    plt.savefig(f"{plot_path}-dist.png")
    plt.clf()

    sns.displot(values, binwidth=1000)
    plt.savefig(f"{plot_path}-dis.png")
    plt.clf()

    for binwidth in [0.01, 0.1]:
        for stat in ['count', 'frequency', 'density']:
            sns.histplot(values, binwidth=binwidth, stat=stat, kde=True)
            plt.savefig(f"{plot_path}-hist-{stat}-{binwidth}.png")
            plt.clf()

    sns.violinplot(values)
    plt.savefig(f"{plot_path}-violin.png")
    plt.clf()

    # plt.show()


def plot_marks_densities(data, output_path):
    output_path = Path(output_path).resolve()
    plot_path = output_path / f"marks-densities"

    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    # Create the data
    collections = list(data)
    collections.sort(reverse=True, key=lambda c: data[c][len(data[c]) - 1] / len(data[c]))
    df = pd.DataFrame(dict(collection=collections))

    # Initialize the FacetGrid object
    # pal = sns.cubehelix_palette(len(df.index), rot=-.25, light=.7)
    pal = sns.color_palette("crest", len(df.index))
    g = sns.FacetGrid(df, row="collection", hue="collection", sharex=False, sharey=False, height=1.75, aspect=10,
                      palette=pal)

    # Draw the densities in a few steps
    bw_adjust = 0.25

    def plot_data(collection, color, label):
        sns.kdeplot(data[label], bw_adjust=bw_adjust, clip_on=False, fill=True, alpha=1, linewidth=1.5, color=color)
        sns.kdeplot(data[label], clip_on=False, color="w", lw=2, bw_adjust=bw_adjust)

    g.map(plot_data, "collection")

    # Passing color=None to refline() uses the hue mapping
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def set_label(x, color, label):
        ax = plt.gca()
        # ax.text(0, .2, label, fontweight="bold", color=color, ha="left", va="center", transform=ax.transAxes)
        ax.text(0, .2, label.upper(), fontsize="x-large", fontweight="bold", color='black', ha="left", va="center",
                transform=ax.transAxes)

    g.map(set_label, "collection")

    # Set the subplots to overlap
    g.figure.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[], xticks=[], ylabel="", xlabel="")
    g.despine(bottom=True, left=True)

    plt.savefig(f"{plot_path}-kde-{bw_adjust}.png")
    plt.savefig(f"{plot_path}-kde-{bw_adjust}.svg")


if __name__ == "__main__":
    main()
