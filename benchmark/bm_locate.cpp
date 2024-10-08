//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
//

#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");
DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

void SetupDefaultCounters(benchmark::State &t_state) {
  t_state.counters["Collection_Size(bytes)"] = 0;
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
  t_state.counters["Occurrences"] = 0;
  t_state.counters["Time_x_Occurrence"] = 0;
}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &_state) {
  for (auto _ : _state) {
    std::vector<int> empty_vector(1000000, 0);
  }

  SetupDefaultCounters(_state);
}
BENCHMARK(BM_WarmUp);

auto BM_QueryLocate = [](benchmark::State &t_state, const auto &t_idx, const auto &t_patterns, auto t_seq_size) {
  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    total_occs = 0;
    for (const auto &pattern: t_patterns) {
      auto occs = t_idx.first->Locate(pattern);
      total_occs += occs.size();
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_seq_size;
  t_state.counters["Size(bytes)"] = t_idx.second;
  t_state.counters["Bits_x_Symbol"] = t_idx.second * 8.0 / t_seq_size;
  t_state.counters["Patterns"] = t_patterns.size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  t_state.counters["Occurrences"] = total_occs;
  t_state.counters["Time_x_Occurrence"] = benchmark::Counter(
      total_occs, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

auto BM_PrintQueryLocate = [](
    benchmark::State &t_state,
    const auto &t_idx_name,
    const auto &t_idx,
    const auto &t_patterns,
    auto t_seq_size) {
  std::string idx_name = t_idx_name;
  replace(idx_name.begin(), idx_name.end(), '/', '_');
  std::string output_filename = "result-" + idx_name + ".txt";

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    std::ofstream out(output_filename);
    total_occs = 0;
    for (const auto &pattern: t_patterns) {
      out << pattern << std::endl;
      auto occs = t_idx.first->Locate(pattern);
      total_occs += occs.size();

      sort(occs.begin(), occs.end());
      for (const auto &item  : occs) {
        out << "  " << item << std::endl;
      }
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_seq_size;
  t_state.counters["Size(bytes)"] = t_idx.second;
  t_state.counters["Bits_x_Symbol"] = t_idx.second * 8.0 / t_seq_size;
  t_state.counters["Patterns"] = t_patterns.size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  t_state.counters["Occurrences"] = total_occs;
  t_state.counters["Time_x_Occurrence"] = benchmark::Counter(
      total_occs, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

int main(int argc, char *argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  std::vector<std::string> patterns;
  {
    std::ifstream pattern_file(FLAGS_patterns.c_str(), std::ios_base::binary);
    if (!pattern_file) {
      std::cerr << "ERROR: Failed to open patterns file! (" << FLAGS_patterns << ")" << std::endl;
      return 3;
    }

    std::string buf;
    while (std::getline(pattern_file, buf)) {
      if (buf.empty())
        continue;

      patterns.emplace_back(buf);
    }
    pattern_file.close();
  }

  // Indexes
  sdsl::cache_config config(true, FLAGS_data_dir, FLAGS_data_name);

  Factory factory(config);

  std::vector<std::pair<const char *, Factory::Config>> index_configs = {
      {"R-Index", Factory::Config{Factory::IndexEnum::RIndex}},

      {"R-Index-WS/4", Factory::Config{Factory::IndexEnum::RIndexSampled, 4}},
      {"R-Index-WS/8", Factory::Config{Factory::IndexEnum::RIndexSampled, 8}},
      {"R-Index-WS/16", Factory::Config{Factory::IndexEnum::RIndexSampled, 16}},
      {"R-Index-WS/32", Factory::Config{Factory::IndexEnum::RIndexSampled, 32}},
      {"R-Index-WS/64", Factory::Config{Factory::IndexEnum::RIndexSampled, 64}},

      {"R-Index-WS-TM/4", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedMarks, 4}},
      {"R-Index-WS-TM/8", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedMarks, 8}},
      {"R-Index-WS-TM/16", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedMarks, 16}},
      {"R-Index-WS-TM/32", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedMarks, 32}},
      {"R-Index-WS-TM/64", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedMarks, 64}},

      {"R-Index-WS-TA/4", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 4}},
      {"R-Index-WS-TA/8", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 8}},
      {"R-Index-WS-TA/16", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 16}},
      {"R-Index-WS-TA/32", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 32}},
      {"R-Index-WS-TA/64", Factory::Config{Factory::IndexEnum::RIndexSampledWithTrustedAreas, 64}},
  };

  std::string print_bm_prefix = "Print-";
  for (const auto &idx_config: index_configs) {
    auto index = factory.make(idx_config.second);

    benchmark::RegisterBenchmark(idx_config.first, BM_QueryLocate, index, patterns, factory.SequenceSize());

    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.first;
      benchmark::RegisterBenchmark(
          print_bm_name.c_str(), BM_PrintQueryLocate, idx_config.first, index, patterns, factory.SequenceSize());
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
