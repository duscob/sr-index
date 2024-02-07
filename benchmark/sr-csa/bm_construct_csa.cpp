//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/10/21.
//

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/memory_management.hpp>

#include "sr-index/construct.h"
#include "sr-index/r_csa.h"
#include "sr-index/sr_csa.h"
#include "sr-index/sr_csa_psi.h"

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

template <typename TIndex>
void BM_ConstructRCSA(benchmark::State &t_state, sri::Config t_config, const std::string &t_data_path) {
  TIndex index;

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-" + idx_name + ".json");
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
}

auto BM_ConstructRCSAWithBWTRun = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructRCSA<sri::RCSAWithBWTRun<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructRCSAWithPsiRun = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructRCSA<sri::RCSAWithPsiRun<>>(t_state, t_config, t_data_path);
};

template<typename TSrIndex>
void BM_ConstructSrIndex(benchmark::State &t_state, sri::Config t_config, const std::string &t_data_path) {
  std::size_t subsample_rate = t_state.range(0); // Subsampling rate

  TSrIndex index(subsample_rate);

  for (auto _ : t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  auto bm_name = t_state.name();
  std::string idx_name = bm_name.substr(bm_name.find('/') + 1);
  {
    std::ofstream ofs("construction-" + idx_name + "-" + std::to_string(subsample_rate) + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-" + idx_name + "-" + std::to_string(subsample_rate) + ".json");
    sdsl::memory_monitor::write_memory_log<sdsl::JSON_FORMAT>(ofs);
    ofs.close();
  }

  setupCommonCounters(t_state);
  t_state.counters["s"] = subsample_rate;
  auto prefix = std::to_string(subsample_rate) + "_";
  {
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(prefix + sri::conf::KEY_BWT_RUN_FIRST_IDX, t_config));
    t_state.counters["r'"] = buf.size();
  }
}

auto BM_ConstructSrCSA = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSA<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSASlim = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSASlim<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSAWithPsiRuns = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSAWithPsiRun<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSAValidMark = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSAValidMark<sri::SrCSA<>>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSAValidArea = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSAValidArea<sri::SrCSA<>>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSRCSAValidMark = [](benchmark::State &t_state, sri::Config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SRCSAValidMark<>>(t_state, t_config, t_data_path);
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

  benchmark::RegisterBenchmark("R-CSA-BWT-Runs", BM_ConstructRCSAWithBWTRun, config, data_path);

  benchmark::RegisterBenchmark("R-CSA-Psi-Runs", BM_ConstructRCSAWithPsiRun, config, data_path);

  benchmark::RegisterBenchmark("SR-CSA", BM_ConstructSrCSA, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-CSA-Slim", BM_ConstructSrCSASlim, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-CSA-Psi-Runs", BM_ConstructSrCSAWithPsiRuns, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-CSA-VM", BM_ConstructSrCSAValidMark, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-CSA-Psi-Runs-VM", BM_ConstructSRCSAValidMark, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::RegisterBenchmark("SR-CSA-VA", BM_ConstructSrCSAValidArea, config, data_path)
      ->RangeMultiplier(2)
      ->Range(FLAGS_min_s, FLAGS_max_s);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
