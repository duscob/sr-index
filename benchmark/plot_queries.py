import sys
import os

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd


def plot_result(_names, _colors, _cmaps, _markers, _sizes, _legend_colors, _title, _xs, _ys, _items, _output_dir):
    scatters = []
    scatters_names = []
    markers = []
    for item in _items:
        plt.plot(_xs[item], _ys[item], color='black', linewidth=0.8)

        colors = _colors
        n = len(_xs[item])
        if n > 1:
            step = 1 / (n - 1)
            colors = [0.0]
            for j in range(n - 1):
                colors.append(colors[j] + step)

        scatters.append(
            plt.scatter(_xs[item], _ys[item], c=colors, marker=_markers[item], s=_sizes[item] / 1.5, cmap=_cmaps[item],
                        vmin=-1., vmax=1., alpha=0.7))  # , label = _names[i]))
        markers.append(plt.plot([], _markers[item], markersize=9, c=_legend_colors[item], label=_names[item]))
        scatters_names.append(_names[item])

    plt.ylim(ymin=0)
    plt.xlim(xmin=0)

    plt.title(_title)

    # plt.yscale('log')
    plt.ylabel('Time (ms / query)')
    plt.xlabel('Size (bps)')

    plt.tight_layout()

    plt.savefig(_output_dir + '/' + 'COLL-' + _title + '-clean.png')
    # plt.legend(scatters, scatters_names)#, bbox_to_anchor=(1.05, 1), loc=2)
    plt.legend()  # edgecolor = 'inherit')
    plt.savefig(_output_dir + '/' + 'COLL-' + _title + '-legend.png')
    plt.show()


if __name__ == "__main__":
    # Count the arguments
    arguments = len(sys.argv) - 1

    if arguments < 4:
        print("The script must be called with: "
              "file of indexes, "
              "file of collections, "
              "directory with query results by collection"
              "output directory")
        exit()

    filename_indexes = sys.argv[1]
    filename_collections = sys.argv[2]
    collection_dir = sys.argv[3]
    output_dir = sys.argv[4]

    print("Indexes:")
    names = {}
    indexes = pd.read_csv(filename_indexes, delimiter=',', encoding='utf-8')
    for i in range(len(indexes.values)):
        print('\t' + indexes['name'][i] + ' => ' + indexes['label'][i])
        names[indexes['name'][i]] = indexes['label'][i]

    print("Collections:")
    collections = pd.read_csv(filename_collections, delimiter=',', encoding='utf-8')
    for c in range(len(collections.values)):
        collection_type = collections['type'][c]
        collection_name = collections['name'][c]
        print('\t' + collection_type + ' - ' + collection_name)

        filename = collection_name + '.csv'
        infile = open(collection_dir + '/' + collection_type + '/' + filename, 'r')
        tmp_filename = 'tmp.csv'

        tmp_file = open(tmp_filename, 'w')
        ln = 0
        for line in infile:
            if ln > 9:
                tmp_file.write(line)
            ln += 1
        tmp_file.close()

        data = pd.read_csv(tmp_filename, delimiter=',', encoding='utf-8')

        if not os.path.exists(output_dir + '/' + collection_type):
            os.makedirs(output_dir + '/' + collection_type)
        outfile = open(output_dir + '/' + collection_type + '/' + filename, 'w')
        outfile.write('label,name,space,time,patterns\n')

        X = []
        Y = []
        indexes_included = {}
        for k in range(len(data.values)):
            index_name = data['name'][k]
            index_name_split = index_name.split('/')
            base_index_name = index_name_split[0]
            if base_index_name in names.keys():
                idx = 0
                if base_index_name not in indexes_included.keys():
                    indexes_included[base_index_name] = len(indexes_included)
                    X.append([])
                    Y.append([])

                idx = indexes_included[base_index_name]
                X[idx].append(data["Bits_x_Symbol"][k])
                Y[idx].append(data["Time_x_Pattern"][k] * 1000)

                outfile.write(names[base_index_name] + ','
                              + index_name + ','
                              + str(X[idx][-1]) + ','
                              + str(Y[idx][-1]) + ','
                              + str(data["Patterns"][k]) + '\n')

        indexes_info = pd.read_csv(collection_dir + '/' + collection_type + '/' + collection_name + '-idx.csv',
                                   delimiter=',', encoding='utf-8')

        colors = [0.5]
        plot_result(indexes['label'], colors, indexes['cmap'], indexes['marker'], indexes['marker_size'],
                    indexes['legend_color'], collections['label'][c], X, Y, indexes_info['show'].values,
                    output_dir + '/' + collection_type)
