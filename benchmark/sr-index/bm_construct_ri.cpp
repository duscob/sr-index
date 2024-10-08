//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/25/22.
//

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/memory_management.hpp>

#include "sr-index/construct.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_string(sa_algo, "SDSL_SE_SAIS", "Suffix Array Algorithm: SDSL_SE_SAIS, SDSL_LIBDIVSUFSORT, BIG_BWT");
DEFINE_int32(min_s, 4, "Minimum sampling parameter s.");
DEFINE_int32(max_s, 2u << 8u, "Maximum sampling parameter s.");

void setupCommonCounters(benchmark::State &t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
}

// Benchmark Warm-up
static void BM_WarmUp(benchmark::State &_state) {
  for (auto _ : _state) {
    std::vector<int> empty_vector(1000000, 0);
  }

  setupCommonCounters(_state);
}
BENCHMARK(BM_WarmUp);

auto BM_ConstructRIndex = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  sri::RIndex<> index;

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  {
    std::ofstream ofs("construction-r-index.html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-r-index.json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  setupCommonCounters(t_state);
  {
    sdsl::int_vector_buffer<8> buf(sdsl::cache_file_name(sdsl::key_text_trait<8>::KEY_TEXT, t_config));
    t_state.counters["n"] = buf.size();
  }
  {
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(sri::conf::KEY_BWT_RUN_FIRST, t_config));
    t_state.counters["r"] = buf.size();
  }
};

template<typename TSrIndex>
void BM_ConstructSRI(benchmark::State &t_state, sri::Config t_config, const std::string &t_data_path) {
  std::size_t subsample_rate = t_state.range(0); // Subsampling rate

  TSrIndex index(subsample_rate);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  {
    std::ofstream ofs("construction-sr-index-" + std::to_string(subsample_rate) + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-sr-index-" + std::to_string(subsample_rate) + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  setupCommonCounters(t_state);
  t_state.counters["s"] = subsample_rate;
  auto prefix = std::to_string(subsample_rate) + "_";
  {
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(prefix + sri::conf::KEY_BWT_RUN_LAST_IDX, t_config));
    t_state.counters["r'"] = buf.size();
  }
};

auto BM_ConstructSRIndex = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSRI<sri::SrIndex<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSRIValidMark = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSRI<sri::SrIndexValidMark<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSRIValidArea = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSRI<sri::SrIndexValidArea<>>(t_state, t_config, t_data_path);
};

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the sr-csa items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  std::string data_path = FLAGS_data;

  sri::Config config(data_path, std::filesystem::current_path(), sri::toSAAlgo(FLAGS_sa_algo));

  benchmark::RegisterBenchmark("R-Index", BM_ConstructRIndex, config, data_path);

  benchmark::RegisterBenchmark("SR-Index", BM_ConstructSRIndex, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-Index-VM", BM_ConstructSRIValidMark, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-Index-VA", BM_ConstructSRIValidArea, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
