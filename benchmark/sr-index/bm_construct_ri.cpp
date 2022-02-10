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
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

void setupCommonCounters(benchmark::State &t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
}

auto BM_ConstructRIndex = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  sri::RIndex<> index;

  for (auto _: t_state) {
    sdsl::memory_monitor::start();
    sri::construct<8>(index, t_data_path, t_config);
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
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(sri::key_trait<8>::KEY_BWT_RUN_FIRST, t_config));
    t_state.counters["r"] = buf.size();
  }
};

template<typename TSrIndex>
void BM_ConstructSrIndex(benchmark::State &t_state, sdsl::cache_config t_config, const std::string &t_data_path) {
  std::size_t subsample_rate = t_state.range(0); // Subsampling rate

  TSrIndex index(subsample_rate);

  for (auto _: t_state) {
    sdsl::memory_monitor::start();
    sri::construct<8>(index, t_data_path, t_config);
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
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(prefix + sri::key_trait<8>::KEY_BWT_RUN_FIRST_IDX, t_config));
    t_state.counters["r'"] = buf.size();
  }
};

auto BM_ConstructSRIndex = [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
  BM_ConstructSrIndex<sri::SRIndex<>>(t_state, t_config, t_data_path);
};

auto BM_ConstructSRIndexValidMark =
    [](benchmark::State &t_state, sdsl::cache_config t_config, const auto &t_data_path) {
      BM_ConstructSrIndex<sri::SRIndexValidMark<>>(t_state, t_config, t_data_path);
    };

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the sr-csa items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::construct_config::byte_algo_sa = FLAGS_sais
                                         ? sdsl::SE_SAIS
                                         : sdsl::LIBDIVSUFSORT; // or LIBDIVSUFSORT for less space-efficient but faster construction

  std::string data_path = FLAGS_data;

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  benchmark::RegisterBenchmark("Construct-R-Index", BM_ConstructRIndex, config, data_path)->Iterations(1);

  benchmark::RegisterBenchmark("Construct-SR-Index", BM_ConstructSRIndex, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::RegisterBenchmark("Construct-SR-Index-ValidMarks", BM_ConstructSRIndexValidMark, config, data_path)
      ->Iterations(1)
      ->RangeMultiplier(2)
      ->Range(4, 2u << 8u);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
