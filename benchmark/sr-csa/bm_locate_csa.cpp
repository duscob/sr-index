//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
//

#include <iostream>

#include <benchmark/benchmark.h>

#include <gflags/gflags.h>

#include <sdsl/config.hpp>

#include "../bm_locate.h"
#include "../bm_locate_set.h"
#include "factory.h"

DEFINE_string(patterns, "", "Patterns file. (MANDATORY)");
DEFINE_string(data_dir, "./", "Data directory.");
DEFINE_string(data_name, "data", "Data file basename.");

DEFINE_int32(min_s, 4, "Minimum sampling parameter s.");
DEFINE_int32(max_s, 128, "Maximum sampling parameter s.");

DEFINE_bool(report_stats, false, "Report statistics for benchmark (mean, median, ...).");
DEFINE_int32(reps, 10, "Repetitions for the locate query benchmark.");
DEFINE_double(min_time, 0, "Minimum time (seconds) for the locate query micro benchmark.");

DEFINE_bool(print_result, false, "Execute benchmark that print results per index.");

int main(int argc, char* argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_patterns.empty() || FLAGS_data_name.empty() || FLAGS_data_dir.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  // Query patterns
  auto patterns = ReadPatterns(FLAGS_patterns);

  // Indexes
  sri::Config config(FLAGS_data_name, FLAGS_data_dir, sri::SDSL_LIBDIVSUFSORT, true);
  auto factory = std::make_shared<Factory<>>(config);

  std::vector<IndexConfig<Factory<>::Config>> idx_configs = {
    //      {"CSA_RAW", Factory<>::Config{Factory<>::IndexEnum::CSA_RAW}, false},
    {"R-CSA-BWT-RUNS", Factory<>::Config{Factory<>::IndexEnum::R_CSA}, false},
    {"R-CSA-PSI-RUNS", Factory<>::Config{Factory<>::IndexEnum::R_CSA_PSI_RUNS}, false},
    {"SR-CSA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA}, true},
    {"SR-CSA-VM", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VM}, true},
    {"SR-CSA-VA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_VA}, true},
    {"SR-CSA-Slim", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM}, true},
    {"SR-CSA-Slim-VM", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VM}, true},
    {"SR-CSA-Slim-VA", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_SLIM_VA}, true},
    {"SR-CSA-PSI-RUNS", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_PSI_RUNS}, true},
    {"SR-CSA-PSI-RUNS-VM", Factory<>::Config{Factory<>::IndexEnum::SR_CSA_PSI_RUNS_VM}, true},
  };

  LocateBenchmarkConfig locate_bm_config{FLAGS_report_stats, FLAGS_reps, FLAGS_min_time, FLAGS_print_result};

  RegisterAllLocateBenchmarks(factory, idx_configs, patterns, locate_bm_config, FLAGS_min_s, FLAGS_max_s);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
