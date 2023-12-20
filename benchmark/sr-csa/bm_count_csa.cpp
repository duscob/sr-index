//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
//

#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

#include "../bm_locate.h"
#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

auto BM_QueryCount = [](benchmark::State &t_state, const auto &t_idx, const auto &t_patterns, auto t_seq_size) {
  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      auto range = t_idx.idx->Count(pattern.decoded);
      total_occs += range.second - range.first;
    }
  }

  UpdateCounter(t_state, t_seq_size, t_idx.size, t_patterns.size(), total_occs);
};

auto BM_PrintQueryCount = [](
    benchmark::State &t_state, const auto &t_idx_name, const auto &t_idx, const auto &t_patterns, auto t_seq_size) {
  std::string idx_name = t_idx_name;
  replace(idx_name.begin(), idx_name.end(), '/', '_');
  std::string output_filename = "result-count-" + idx_name + ".csv";

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    std::ofstream out(output_filename);
    out << "pattern,count,range_start,range_end" << std::endl;
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      auto range = t_idx.idx->Count(pattern.decoded);
      auto count = range.second - range.first;
      total_occs += count;
      out << "\"" << pattern.encoded << "\"," << count << "," << range.first << "," << range.second << std::endl;
    }
  }

  UpdateCounter(t_state, t_seq_size, t_idx.size, t_patterns.size(), total_occs);
};

int main(int argc, char *argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  auto patterns = ReadPatterns(FLAGS_patterns);

  // Indexes
  const sri::Config config(FLAGS_data_name, FLAGS_data_dir, sri::SDSL_LIBDIVSUFSORT, true);

  Factory<> factory(config);

  std::vector<std::pair<const char *, Factory<>::Config>> index_configs = {
      {"R-CSA", Factory<>::Config{Factory<>::IndexEnum::R_CSA}},
  };

  std::string print_bm_prefix = "Print-";
  for (const auto &idx_config : index_configs) {
    auto index = factory.make(idx_config.second);

    benchmark::RegisterBenchmark(idx_config.first, BM_QueryCount, index, patterns, factory.sizeSequence());

    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.first;
      benchmark::RegisterBenchmark(
          print_bm_name, BM_PrintQueryCount, idx_config.first, index, patterns, factory.sizeSequence());
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
