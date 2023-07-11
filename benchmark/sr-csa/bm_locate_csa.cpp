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
DEFINE_int32(min_s, 4, "Minimum sampling parameter s.");
DEFINE_int32(max_s, 128, "Maximum sampling parameter s.");
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

auto MakeIndex = [](auto t_factory, auto &t_config, auto &t_state) {
  if (t_config.has_sampling) {
    t_config.fac_config.sampling_size = t_state.range(0);
  }

  return t_factory->make(t_config.fac_config);
};

auto BM_QueryLocate = [](benchmark::State &t_state, auto t_factory, auto t_config, const auto &t_patterns) {
  auto idx = MakeIndex(t_factory, t_config, t_state);

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      auto occs = idx.first->Locate(pattern);
      total_occs += occs.size();
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_factory->sizeSequence();
  t_state.counters["Size(bytes)"] = idx.second;
  t_state.counters["Bits_x_Symbol"] = idx.second * 8.0 / t_factory->sizeSequence();
  t_state.counters["Patterns"] = t_patterns.size();
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_patterns.size(), benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  t_state.counters["Occurrences"] = total_occs;
  t_state.counters["Time_x_Occurrence"] = benchmark::Counter(
      total_occs, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

auto BM_PrintQueryLocate = [](benchmark::State &t_state, auto t_factory, auto t_config, const auto &t_patterns) {
  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  replace(idx_name.begin(), idx_name.end(), '/', '_');
  std::string output_filename = "result-" + idx_name + ".txt";

  auto idx = MakeIndex(t_factory, t_config, t_state);

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    std::ofstream out(output_filename);
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      out << pattern << std::endl;
      auto occs = idx.first->Locate(pattern);
      total_occs += occs.size();

      sort(occs.begin(), occs.end());
      for (const auto &item : occs) {
        out << "  " << item << std::endl;
      }
    }
  }

  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_factory->sizeSequence();
  t_state.counters["Size(bytes)"] = idx.second;
  t_state.counters["Bits_x_Symbol"] = idx.second * 8.0 / t_factory->sizeSequence();
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

  auto factory = std::make_shared<Factory<>>(config);

  struct BenchmarkConfig {
    std::string name;
    Factory<>::Config fac_config;
    bool has_sampling = false;
  };
  std::vector<BenchmarkConfig> bm_configs = {
//      {"CSA_RAW", Factory<>::Config{Factory<>::IndexEnum::CSA_RAW}, false},
      {"R-CSA", Factory<>::Config{Factory<>::IndexEnum::R_CSA}, false},
      {"SR-CSA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA}, true},
      {"SR-CSA-VM", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM}, true},
      {"SR-CSA-VA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA}, true},
      {"SR-CSA-Slim", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM}, true},
      {"SR-CSA-Slim-VM", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM}, true},
      {"SR-CSA-Slim-VA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA}, true},
  };

  std::string print_bm_prefix = "Print/";
  for (const auto &bm_config : bm_configs) {
    try {
      auto bm = benchmark::RegisterBenchmark(bm_config.name, BM_QueryLocate, factory, bm_config, patterns);
      if (bm_config.has_sampling) {
        bm->RangeMultiplier(2)->Range(FLAGS_min_s, FLAGS_max_s);
      }

      if (FLAGS_print_result) {
        auto print_bm_name = print_bm_prefix + bm_config.name;
        auto bm_print = benchmark::RegisterBenchmark(print_bm_name, BM_PrintQueryLocate, factory, bm_config, patterns);
        if (bm_config.has_sampling) {
          bm_print->RangeMultiplier(2)->Range(FLAGS_min_s, FLAGS_max_s);
        }
      }
    }
    catch (const std::exception &error) {
      std::cerr << error.what() << std::endl;
    }
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
