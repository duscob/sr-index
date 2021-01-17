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
  t_state.counters["Size(bytes)"] = 0;
  t_state.counters["Bits_x_Symbol"] = 0;
  t_state.counters["Patterns"] = 0;
  t_state.counters["Time_x_Pattern"] = 0;
}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &_state) {
  for (auto _ : _state) {
    std::string empty_string;
  }

  SetupDefaultCounters(_state);
}
BENCHMARK(BM_WarmUp);

// Benchmark Queries on Document Frequency Index
auto BM_QueryLocate =
    [](benchmark::State &_state, const auto &_idx, const auto &_patterns, auto _seq_size) {
      for (auto _ : _state) {
        for (const auto &pattern: _patterns) {
          auto occs = _idx.first->Locate(pattern);
        }
      }

      SetupDefaultCounters(_state);
      _state.counters["Size(bytes)"] = _idx.second;
      _state.counters["Bits_x_Symbol"] = _idx.second * 8.0 / _seq_size;
      _state.counters["Patterns"] = _patterns.size();
      _state.counters["Time_x_Pattern"] = benchmark::Counter(
          _patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
    };

auto BM_PrintQueryLocate =
    [](benchmark::State &_state, const auto &_idx_name, const auto &_idx, const auto &_patterns, auto _seq_size) {
      std::string idx_name = _idx_name;
      replace(idx_name.begin(), idx_name.end(), '/', '_');
      std::string output_filename = "result-" + idx_name + ".txt";

      for (auto _ : _state) {
        std::ofstream out(output_filename);
        for (const auto &pattern: _patterns) {
          out << pattern << std::endl;
          auto occs = _idx.first->Locate(pattern);

          for (const auto &item  : occs) {
            out << "  " << item << std::endl;
          }
        }
      }

      SetupDefaultCounters(_state);
      _state.counters["Size(bytes)"] = _idx.second;
      _state.counters["Bits_x_Symbol"] = _idx.second * 8.0 / _seq_size;
      _state.counters["Patterns"] = _patterns.size();
      _state.counters["Time_x_Pattern"] = benchmark::Counter(
          _patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
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
      std::cerr << "ERROR: Failed to open patterns file!" << std::endl;
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