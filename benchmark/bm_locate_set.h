//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 7/13/2023.
//

#ifndef SRI_BENCHMARK_BM_LOCATE_SET_H_
#define SRI_BENCHMARK_BM_LOCATE_SET_H_

#include <string>

#include "bm_locate.h"

template<typename TFactoryConfig>
struct IndexConfig {
  std::string name;
  TFactoryConfig factory_config;
  bool has_sampling = false;
};

auto RegisterAllLocateBenchmarks = [](
    auto t_factory, const auto &t_idx_configs, const auto &t_patterns, auto t_bm_config, auto t_min_s, auto t_max_s
) {
  auto n = t_factory->sizeSequence();

  for (auto idx_config : t_idx_configs) {
    // Index builder
    auto make_index = [t_factory, idx_config](const auto &tt_state) {
      auto fac_config = idx_config.factory_config;
      if (idx_config.has_sampling) {
        fac_config.sampling_size = tt_state.range(0);
      }

      auto idx = t_factory->make(fac_config);

      auto locate = [idx](const auto &ttt_pattern) {
        return idx.idx->Locate(ttt_pattern);
      };

      return std::make_pair(locate, idx.size);
    };

    // Benchmarks
    auto bms =
        RegisterLocateBenchmarks(idx_config.name, make_index, t_patterns, n, t_bm_config, idx_config.has_sampling);
    // Index with sampling
    if (idx_config.has_sampling) {
      for (auto &bm : bms) {
        bm.second->RangeMultiplier(2)->Range(t_min_s, t_max_s);
      }
    }
  }
};

#endif //SRI_BENCHMARK_BM_LOCATE_SET_H_
