//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/12/2023.
//

#ifndef SRI_BENCHMARK_BM_LOCATE_H_
#define SRI_BENCHMARK_BM_LOCATE_H_

#include <fstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <vector>

#include <benchmark/benchmark.h>

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

auto UpdateCounter = [](benchmark::State &t_state, auto t_n, auto t_index_size, auto t_n_patterns, auto t_n_occs) {
  SetupDefaultCounters(t_state);
  t_state.counters["Collection_Size(bytes)"] = t_n;
  t_state.counters["Size(bytes)"] = t_index_size;
  t_state.counters["Bits_x_Symbol"] = t_index_size * 8.0 / t_n;
  t_state.counters["Patterns"] = t_n_patterns;
  t_state.counters["Time_x_Pattern"] = benchmark::Counter(
      t_n_patterns, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  t_state.counters["Occurrences"] = t_n_occs;
  t_state.counters["Time_x_Occurrence"] = benchmark::Counter(
      t_n_occs, benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
};

auto BM_MacroLocate = [](benchmark::State &t_state, auto t_make_index, const auto &t_patterns, auto t_n) {
  auto [locate, index_size] = t_make_index(t_state);

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      auto occs = locate(pattern);
      total_occs += occs.size();
    }
  }

  UpdateCounter(t_state, t_n, index_size, t_patterns.size(), total_occs);
};

auto BM_MicroLocate = [](benchmark::State &t_state, auto t_make_index, const auto &t_patterns, auto t_i, auto t_n) {
  auto [locate, index_size] = t_make_index(t_state);

  const auto &pattern = t_patterns[*t_i];
  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    auto occs = locate(pattern);
    total_occs = occs.size();
  }

  if (++*t_i == t_patterns.size()) {
    *t_i = 0;
  }

  UpdateCounter(t_state, t_n, index_size, 1, total_occs);
};

auto BM_PrintLocate = [](
    benchmark::State &t_state, auto t_make_index, const auto &t_patterns, auto t_n, bool t_is_sampled_index
) {
  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  replace(idx_name.begin(), idx_name.end(), '/', '_');
  if (t_is_sampled_index) {
    idx_name += "-" + std::to_string(t_state.range(0));
  }
  std::string output_filename = "result-" + idx_name + ".txt";

  auto [locate, index_size] = t_make_index(t_state);

  std::size_t total_occs = 0;

  for (auto _ : t_state) {
    std::ofstream out(output_filename);
    total_occs = 0;
    for (const auto &pattern : t_patterns) {
      out << pattern << std::endl;
      auto occs = locate(pattern);
      total_occs += occs.size();

      sort(occs.begin(), occs.end());
      for (const auto &item : occs) {
        out << "  " << item << std::endl;
      }
    }
  }

  UpdateCounter(t_state, t_n, index_size, t_patterns.size(), total_occs);
};

enum KeyLocateBenchmark {
  kMacro,
  kMicro,
  kPrint
};

using LocateBenchmarks = std::map<KeyLocateBenchmark, benchmark::internal::Benchmark *>;

struct LocateBenchmarkConfig {
  bool report_stats = false;
  int reps = 10;
  double min_time = 0.0;
  bool print_results = false;
};

auto RegisterLocateBenchmarks = [](
    const auto &t_name,
    auto t_make_index,
    const auto &t_patterns,
    auto t_n,
    const LocateBenchmarkConfig &t_bm_config,
    bool t_is_sampled_index
) {
  LocateBenchmarks bms;

  bms[kMacro] = benchmark::RegisterBenchmark(t_name, BM_MacroLocate, t_make_index, t_patterns, t_n);

  if (t_bm_config.report_stats) {
    auto statistics_min = [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
    };
    auto statistics_max = [](const std::vector<double> &v) -> double {
      return *(std::max_element(std::begin(v), std::end(v)));
    };

    bms[kMacro]->Name(t_name + "/macro")
        ->Repetitions(t_bm_config.reps)
        ->ComputeStatistics("min", statistics_min)
        ->ComputeStatistics("max", statistics_max)
        ->ReportAggregatesOnly();

    bms[kMicro] = benchmark::RegisterBenchmark(
        t_name + "/micro", BM_MicroLocate, t_make_index, t_patterns, std::make_shared<int>(0), t_n)
        ->Repetitions(t_patterns.size())
        ->ComputeStatistics("min", statistics_min)
        ->ComputeStatistics("max", statistics_max)
        ->ComputeStatistics("total", [](const std::vector<double> &v) -> double {
          return std::accumulate(std::begin(v), std::end(v), double(0));
        })
        ->ReportAggregatesOnly();
    if (t_bm_config.min_time > 0) {
      bms[kMicro]->MinTime(t_bm_config.min_time);
    } else {
      bms[kMicro]->Iterations(t_bm_config.reps);
    }
  }

  if (t_bm_config.print_results) {
    bms[kPrint] = benchmark::RegisterBenchmark(
        "Print/" + t_name, BM_PrintLocate, t_make_index, t_patterns, t_n, t_is_sampled_index)
        ->Iterations(1);
  }

  return bms;
};

std::vector<std::string> ReadPatterns(const std::string &t_pattern_path) {
  std::vector<std::string> patterns;

  std::ifstream pattern_file(t_pattern_path, std::ios_base::binary);
  if (!pattern_file) {
    std::cerr << "ERROR: Failed to open patterns file! (" << t_pattern_path << ")" << std::endl;
    exit(3);
  }

  std::string buffer;
  while (std::getline(pattern_file, buffer)) {
    if (buffer.empty())
      continue;

    patterns.emplace_back(buffer);
  }
  pattern_file.close();

  return patterns;
}

#endif //SRI_BENCHMARK_BM_LOCATE_H_
