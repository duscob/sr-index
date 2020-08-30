#!/usr/bin/env bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

dir_colls=$1
dir_tools=${2:-"$SCRIPT_DIR/../build"}

for coll in "$dir_colls"/*; do
  coll_name=$(basename "$coll")
  echo "Collection $coll_name"

  mkdir -p "$coll_name"
  cd "$coll_name" || exit

  "$dir_tools"/bm_build_ri_items --data "$coll"/data --benchmark_counters_tabular=true --benchmark_out_format=csv --benchmark_out="$coll_name"-build.csv

  cd ..
done
