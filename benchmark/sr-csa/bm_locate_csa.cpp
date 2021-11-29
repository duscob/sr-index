//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
//

#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

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
  for (auto _: _state) {
    std::vector<int> empty_vector(1000000, 0);
  }

  SetupDefaultCounters(_state);
}
BENCHMARK(BM_WarmUp);

auto BM_QueryLocate = [](benchmark::State &t_state, const auto &t_idx, const auto &t_patterns, auto t_seq_size) {
  std::size_t total_occs = 0;

  for (auto _: t_state) {
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
    benchmark::State &t_state, const auto &t_idx_name, const auto &t_idx, const auto &t_patterns, auto t_seq_size) {
  std::string idx_name = t_idx_name;
  replace(idx_name.begin(), idx_name.end(), '/', '_');
  std::string output_filename = "result-" + idx_name + ".txt";

  std::size_t total_occs = 0;

  for (auto _: t_state) {
    std::ofstream out(output_filename);
    total_occs = 0;
    for (const auto &pattern: t_patterns) {
      out << pattern << std::endl;
      auto occs = t_idx.first->Locate(pattern);
      total_occs += occs.size();

      sort(occs.begin(), occs.end());
      for (const auto &item: occs) {
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

  Factory<> factory(config);

  std::vector<std::pair<const char *, Factory<>::Config>> index_configs = {
      {"CSA_RAW", Factory<>::Config{Factory<>::IndexEnum::CSA_RAW}},

      {"R-CSA", Factory<>::Config{Factory<>::IndexEnum::R_CSA}},

      {"SR-CSA/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA, 4}},
      {"SR-CSA/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA, 8}},
      {"SR-CSA/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA, 16}},
      {"SR-CSA/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA, 32}},
      {"SR-CSA/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA, 64}},

      {"SR-CSA-VM/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM, 4}},
      {"SR-CSA-VM/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM, 8}},
      {"SR-CSA-VM/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM, 16}},
      {"SR-CSA-VM/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM, 32}},
      {"SR-CSA-VM/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM, 64}},

      {"SR-CSA-VA/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 4}},
      {"SR-CSA-VA/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 8}},
      {"SR-CSA-VA/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 16}},
      {"SR-CSA-VA/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 32}},
      {"SR-CSA-VA/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 64}},
//      {"SR-CSA-VA/128", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 128}},
//      {"SR-CSA-VA/256", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA, 256}},

      {"SR-CSA-SLIM/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM, 4}},
      {"SR-CSA-SLIM/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM, 8}},
      {"SR-CSA-SLIM/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM, 16}},
      {"SR-CSA-SLIM/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM, 32}},
      {"SR-CSA-SLIM/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM, 64}},

      {"SR-CSA-SLIM-VM/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM, 4}},
      {"SR-CSA-SLIM-VM/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM, 8}},
      {"SR-CSA-SLIM-VM/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM, 16}},
      {"SR-CSA-SLIM-VM/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM, 32}},
      {"SR-CSA-SLIM-VM/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM, 64}},

      {"SR-CSA-SLIM-VA/4", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 4}},
      {"SR-CSA-SLIM-VA/8", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 8}},
      {"SR-CSA-SLIM-VA/16", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 16}},
      {"SR-CSA-SLIM-VA/32", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 32}},
      {"SR-CSA-SLIM-VA/64", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 64}},
//      {"SR-CSA-SLIM-VA/128", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 128}},
//      {"SR-CSA-SLIM-VA/256", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA, 256}},
  };

  std::string print_bm_prefix = "Print-";
  for (const auto &idx_config: index_configs) {
    auto index = factory.make(idx_config.second);

    benchmark::RegisterBenchmark(idx_config.first, BM_QueryLocate, index, patterns, factory.sizeSequence());

    if (FLAGS_print_result) {
      auto print_bm_name = print_bm_prefix + idx_config.first;
      benchmark::RegisterBenchmark(
          print_bm_name.c_str(), BM_PrintQueryLocate, idx_config.first, index, patterns, factory.sizeSequence());
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
