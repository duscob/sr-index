#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

dir_idxs=$1
dir_colls=$2
benchmark=${3:-"$SCRIPT_DIR/../build/bm_locate"}

for coll in "$dir_idxs"/*; do
  coll_name=$(basename "$coll")
  echo "Collection $coll_name"

  mkdir -p "$coll_name"
  cd "$coll_name" || exit

  $benchmark --data_dir "$dir_idxs"/"$coll_name" --patterns "$dir_colls"/"$coll_name"/patterns --benchmark_counters_tabular=true --benchmark_out_format=csv --benchmark_out="$coll_name".csv --print_result

  cd ..
done
