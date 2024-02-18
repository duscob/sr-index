import argparse
from pathlib import Path
import subprocess
import struct

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

    collection_path = Path(args.collection)
    if args.group:
        # Processing collections group
        for child in collection_path.iterdir():
            if child.is_dir():
                process_collection(child, args.output_path, args.cmd_path)
    else:
        # Processing single collection
        process_collection(collection_path, args.output_path, args.cmd_path)


def esc(code):
    return f'\033[{code}m'


def process_collection(collection_path, output_path, bm_cmd_path):
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



if __name__ == "__main__":
    main()
