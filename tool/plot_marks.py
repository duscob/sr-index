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

    data = {"inter": {}, "intra": {}}
    collection_path = Path(args.collection)
    if args.group:
        # Processing collections group
        for child in collection_path.iterdir():
            if child.is_dir():
                process_collection(child, args.output_path, args.cmd_path, data)
    else:
        # Processing single collection
        process_collection(collection_path, args.output_path, args.cmd_path, data)

    plot_marks_densities_grid(data["inter"], args.output_path)

    for key, value in data["intra"].items():
        plot_marks_densities(key, value, args.output_path)
        plot_marks_frequencies(key, value, args.output_path)


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
    # output_path.mkdir(parents=True, exist_ok=True)

    data["inter"][collection_name] = process_data(bm_cmd_path, collection_path, "bwt_run_first_text_pos")
    data["intra"][collection_name] = {}
    data["intra"][collection_name]["0"] = data["inter"][collection_name]

    for sampling_factor in ["4", "8", "16", "32", "64"]:
        data["intra"][collection_name][sampling_factor] = process_data(
            bm_cmd_path, collection_path, f"{sampling_factor}_bwt_run_first_text_pos")

    # plot_marks_density(values, collection_name, output_path)


def process_data(bm_cmd_path, collection_path, data_key):
    # data_path = collection_path / "sri" / "bwt_run_first_text_pos_data.sdsl.vec.bin"
    data_path = collection_path / "sri" / f"{data_key}_data.sdsl.vec.bin"
    if not data_path.exists():
        print(f"{esc('38;5;69')}Create file for marks vector ({data_path}) {esc(0)}")

        cmd = str(bm_cmd_path / "int_vector_to_vector")
        cmd += " --data=./data"
        cmd += f" --key={data_key}"
        cmd += " 2>int_vector_to_vector-error.txt"

        subprocess.run(cmd, shell=True, cwd=collection_path / "sri")
    values = read_data(data_path)
    # values = np.array(read_data(data_path)).astype(np.float64)
    values.sort()
    # values *= 100.0 / values[len(values) - 1]

    # data["inter"][collection_name] = values
    return values


def read_data(data_path):
    values = []
    with open(data_path, 'rb') as fin:
        for chunk in iter(lambda: fin.read(8), b''):
            values.append(struct.unpack('Q', chunk)[0])

    return values


def plot_marks_density(values, collection_name, output_path):
    output_path.mkdir(parents=True, exist_ok=True)
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


def plot_marks_densities_grid(data, output_path):
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
        # sns.histplot(data[label], clip_on=False, stat="frequency", bins=int(len(data[label]) / 1000), facecolor=color,
        #              kde=True, kde_kws={"bw_adjust": bw_adjust}, color="red")
        # sns.displot(data[label], kde=True, color='red',
        #             line_kws={'lw': 3}, facecolor=color, edgecolor='black')
        # print(len(data[label]))
        # sns.displot(data[label], kind="kde", fill=True)

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
    # plt.savefig(f"{plot_path}-kde-{bw_adjust}.svg")
    plt.clf()


def plot_marks_densities(name, data, output_path):
    output_path = Path(output_path).resolve()
    plot_path = output_path / f"marks-densities"

    style = {
        "axes.facecolor": (0, 0, 0, 0),
        "axes.titlesize": "xx-large",
        "axes.labelsize": "x-large",
        "xtick.labelsize": "large",
        "ytick.labelsize": "large"
    }
    sns.set_theme(style="white", rc=style)

    # Create the data
    collections = list(data)
    collections.sort(reverse=True, key=lambda c: data[c][len(data[c]) - 1] / len(data[c]))
    df = pd.DataFrame(dict(collection=collections))

    pal = sns.color_palette("crest", len(df.index))
    # pal = sns.color_palette("crest", len(df.index))

    bw_adjust = 0.25

    fig = sns.displot(data, kind="kde", bw_adjust=bw_adjust, clip_on=False, fill=True, aspect=4, palette="crest",
                      legend=False)

    # sns.kdeplot(data, bw_adjust=bw_adjust, clip_on=False, fill=True, alpha=0.5, linewidth=0.5, aspect=4, legend=False).set_title(name)

    fig.set(title=f"{name}")
    # fig.set(xlim=0)

    plt.tight_layout(pad=0.5)
    plt.savefig(f"{plot_path}-kde-{bw_adjust}-{name}.png")
    # plt.savefig(f"{plot_path}-kde-{bw_adjust}-{name}.svg")
    plt.clf()


def plot_marks_frequencies(name, data, output_path):
    output_path = Path(output_path).resolve()
    plot_path = output_path / f"marks-densities"

    style = {
        "axes.facecolor": (0, 0, 0, 0),
        "axes.titlesize": "xx-large",
        "axes.labelsize": "x-large",
        "xtick.labelsize": "large",
        "ytick.labelsize": "large"
    }
    sns.set_theme(style="white", rc=style)

    # Create the data
    collections = list(data)
    collections.sort(reverse=True, key=lambda c: data[c][len(data[c]) - 1] / len(data[c]))
    # df = pd.DataFrame(dict(collection=collections))
    #
    # pal = sns.color_palette("crest", len(df.index))

    # bw_adjust = 0.25
    # (sns
    #  .displot(data, stat="frequency", kind="hist", kde=True, kde_kws={"bw_adjust":bw_adjust}, clip_on=False, fill=True, aspect=4, palette="crest", legend=False)
    #  .set(title=f"{name}")
    #  )
    fig = sns.displot(data, stat="frequency", kind="hist", element="poly", clip_on=False, fill=True, aspect=4,
                      palette="crest", legend=False)

    fig.set(title=f"{name}")
    fig.set(xlim=0)

    plt.tight_layout(pad=0.5)

    plt.savefig(f"{plot_path}-freq-{name}.png")
    # plt.savefig(f"{plot_path}-freq-{name}.svg")
    plt.clf()

    ranges = [(0, data["0"][len(data["0"]) - 1] / 4),
              (data["0"][len(data["0"]) - 1] / 4, data["0"][len(data["0"]) - 1] / 2),
              (data["0"][len(data["0"]) - 1] / 2, data["0"][len(data["0"]) - 1] * 3 / 4),
              (data["0"][len(data["0"]) - 1] * 3 / 4, data["0"][len(data["0"]) - 1])
              ]
    i = 1
    for l, h in ranges:
        # data_filtered = {}
        # for key, value in data.items():
        #     data_filtered[key] = [v for v in value if l <= v <= h]

        sns.set_theme(style="white", font_scale=2.5, rc=style)
        fig = sns.displot(data, stat="frequency", kind="hist", element="poly", clip_on=False, fill=True,
                          height=40, aspect=0.67, palette="crest", legend=False)

        fig.set(title=f"{name}")
        fig.set(xlim=(l, h))

        plt.tight_layout(pad=0.5)
        plt.savefig(f"{plot_path}-freq-{name}-part-{i}.png")
        plt.clf()

        i = i + 1


if __name__ == "__main__":
    main()
