//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_FACTORY_H_
#define SRI_BENCHMARK_SR_CSA_FACTORY_H_

#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "sr-index/config.h"
#include "sr-index/r_csa.h"
#include "sr-index/r_csa_bwt.h"
#include "sr-index/sr_csa.h"
#include "sr-index/sr_csa_bwt.h"

#include "../bm_base.h"

class Factory {
 public:
  enum class IndexEnum {
    R_CSA = 0,
    SR_CSA,
    SR_CSA_VM,
    SR_CSA_VA,
    R_CSA_BWT,
    CSA_RAW,
    SR_CSA_BWT,
    SR_CSA_BWT_VM,
    SR_CSA_BWT_VA,
    SR_CSA_BWT_SLIM,
    SR_CSA_BWT_VM_SLIM,
    SR_CSA_BWT_VA_SLIM,
  };

  struct Config {
    IndexEnum index_type;
    std::size_t sampling_size;

    bool operator<(const Config& t_c) const {
      return index_type < t_c.index_type || (index_type == t_c.index_type && sampling_size < t_c.sampling_size);
    }
  };

  explicit Factory(sri::Config t_config) : config_{std::move(t_config)} {
    using namespace sri::conf;
    sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(config_.keys[kBWT][kBase], config_));
    n_ = buf.size();
  }

  auto sizeSequence() const {
    return n_;
  }

  struct Index {
    std::shared_ptr<sri::LocateIndex<>> idx;
    std::size_t size = 0;
  };

  // Count-only index: loads strictly fewer components than the locate variants. Useful for
  // measuring the load-time / memory cost of count-only queries against the same corpus.
  struct CountIndex {
    std::shared_ptr<sri::CountIndex<>> idx;
    std::size_t size = 0;
  };

  CountIndex makeCount() {
    if (count_index_.idx) {
      return count_index_;
    }
    auto idx = std::make_shared<sri::RCSACount<ExternalGenericStorage>>(std::ref(storage_));
    idx->load(config_);
    count_index_ = {idx, sdsl::size_in_bytes(*idx)};
    return count_index_;
  }

  Index make(const Config& t_config) {
    auto it = indexes_.find(t_config);
    if (it != indexes_.end()) {
      return it->second;
    }

    Index index;
    switch (t_config.index_type) {
      case IndexEnum::R_CSA: {
        auto idx = std::make_shared<sri::RCSA<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA: {
        auto idx = std::make_shared<sri::SrCSA<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_VM: {
        auto idx = std::make_shared<sri::SrCSAValidMark<sri::SrCSA<ExternalGenericStorage>>>(std::ref(storage_),
                                                                                             t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_VA: {
        auto idx = std::make_shared<sri::SrCSAValidArea<sri::SrCSA<ExternalGenericStorage>>>(std::ref(storage_),
                                                                                             t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::R_CSA_BWT: {
        auto idx = std::make_shared<sri::RCSABWTRun<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::CSA_RAW: {
        auto idx = std::make_shared<sri::CSARaw<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT: {
        auto idx =
            std::make_shared<sri::SrCSABWTRun<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT_VM: {
        auto idx = std::make_shared<sri::SrCSABWTRunValidMark<sri::SrCSABWTRun<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT_VA: {
        auto idx = std::make_shared<sri::SrCSABWTRunValidArea<sri::SrCSABWTRun<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT_SLIM: {
        auto idx =
            std::make_shared<sri::SrCSABWTRunSlim<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT_VM_SLIM: {
        auto idx = std::make_shared<sri::SrCSABWTRunValidMark<sri::SrCSABWTRunSlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_BWT_VA_SLIM: {
        auto idx = std::make_shared<sri::SrCSABWTRunValidArea<sri::SrCSABWTRunSlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }
    }

    if (index.idx) {
      indexes_[t_config] = index;
    }

    return index;
  }

 private:
  sri::Config config_;

  std::size_t n_ = 0;

  sri::GenericStorage storage_;

  std::map<Config, Index> indexes_;
  CountIndex count_index_;
};

#endif  // SRI_BENCHMARK_SR_CSA_FACTORY_H_
