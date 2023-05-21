//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/10/21.
//

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/memory_management.hpp>

#include "sr-index/construct.h"
#include "sr-index/csa.h"
#include "sr-index/sr_csa.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
//DEFINE_bool(rebuild, false, "Rebuild all the items.");
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

void setupCommonCounters(benchmark::State &t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
//  t_state.counters["mr'"] = 0;
}

auto BM_ConstructCSA = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  sri::CSA<> index;

  for (auto _: t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  {
    std::ofstream ofs("construction-r-csa.html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-r-csa.json");
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
void BM_ConstructSrIndex(benchmark::State &t_state, sdsl::cache_config t_config, const std::string &t_data_path) {
  std::size_t subsample_rate = t_state.range(0); // Subsampling rate

  TSrIndex index(subsample_rate);

  for (auto _: t_state) {
    sdsl::memory_monitor::start();
    sri::construct(index, t_data_path, t_config);
    sdsl::memory_monitor::stop();
  }

  {
    std::ofstream ofs("construction-sr-csa-" + std::to_string(subsample_rate) + ".html");
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(ofs);
    ofs.close();
  }
  {
    std::ofstream ofs("construction-sr-csa-" + std::to_string(subsample_rate) + ".json");
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
};

auto BM_ConstructSrCSA = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSA<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSASlim = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSASlim<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSAValidMark = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSAValidMark<sri::SrCSA<>>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSrCSAValidArea = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SrCSAValidArea<sri::SrCSA<>>>(t_state, t_config, t_data_path);
};

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the sr-csa items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::construct_config().byte_algo_sa = FLAGS_sais
                                         ? sdsl::SE_SAIS
                                         : sdsl::LIBDIVSUFSORT; // or LIBDIVSUFSORT for less space-efficient but faster construction

  std::string data_path = FLAGS_data;

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  benchmark::RegisterBenchmark("Construct-R-CSA", BM_ConstructCSA, config, data_path)->Iterations(1);

  benchmark::RegisterBenchmark("Construct-SR-CSA", BM_ConstructSrCSA, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::RegisterBenchmark("Construct-SR-CSA-Slim", BM_ConstructSrCSASlim, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::RegisterBenchmark("Construct-SR-CSA-ValidMark", BM_ConstructSrCSAValidMark, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::RegisterBenchmark("Construct-SR-CSA-ValidArea", BM_ConstructSrCSAValidArea, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
