//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_FACTORY_H_
#define SRI_BENCHMARK_SR_CSA_FACTORY_H_

#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "sr-index/r_csa.h"
#include "sr-index/sr_csa.h"
#include "sr-index/sr_csa_psi.h"
#include "config.h"

using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

template<uint8_t t_width = 8>
class Factory {
 public:
  enum class IndexEnum {
    R_CSA = 0,
    CSA_RAW,
    SR_CSA,
    SR_CSA_VM,
    SR_CSA_VA,
    SR_CSA_SLIM,
    SR_CSA_SLIM_VM,
    SR_CSA_SLIM_VA,
    R_CSA_PSI_RUNS,
    SR_CSA_PSI_RUNS,
    SR_CSA_PSI_RUNS_VM,
    SR_CSA_PSI_RUNS_VA
  };

  struct Config {
    IndexEnum index_type;
    std::size_t sampling_size;

    bool operator<(const Config &t_c) const {
      return index_type < t_c.index_type || (index_type == t_c.index_type && sampling_size < t_c.sampling_size);
    }
  };

  explicit Factory(sri::Config t_config) : config_{std::move(t_config)} {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
    n_ = buf.size();
  }

  auto sizeSequence() const { return n_; }

  struct Index {
    std::shared_ptr<sri::LocateIndex> idx;
    std::size_t size = 0;
  };

  Index make(const Config &t_config) {
    auto it = indexes_.find(t_config);
    if (it != indexes_.end()) {
      return it->second;
    }

    Index index;
    switch (t_config.index_type) {
      case IndexEnum::R_CSA: {
        auto idx = std::make_shared<sri::RCSAWithBWTRun<ExternalGenericStorage>>(std::ref(storage_));
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

      case IndexEnum::SR_CSA: {
        auto idx = std::make_shared<sri::SrCSA<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_VM: {
        auto idx = std::make_shared<sri::SrCSAValidMark<sri::SrCSA<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_VA: {
        auto idx = std::make_shared<sri::SrCSAValidArea<sri::SrCSA<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_SLIM: {
        auto idx = std::make_shared<sri::SrCSASlim<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_SLIM_VM: {
        auto idx = std::make_shared<sri::SrCSAValidMark<sri::SrCSASlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_SLIM_VA: {
        auto idx = std::make_shared<sri::SrCSAValidArea<sri::SrCSASlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::R_CSA_PSI_RUNS: {
        auto idx = std::make_shared<sri::RCSAWithPsiRun<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_PSI_RUNS: {
        auto idx = std::make_shared<sri::SrCSAWithPsiRun<ExternalGenericStorage>>(std::ref(storage_),
          t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_PSI_RUNS_VM: {
        auto idx = std::make_shared<sri::SRCSAValidMark<sri::SrCSAWithPsiRun<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_CSA_PSI_RUNS_VA: {
        auto idx = std::make_shared<sri::SRCSAValidArea<sri::SrCSAWithPsiRun<ExternalGenericStorage>>>(
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
};

#endif //SRI_BENCHMARK_SR_CSA_FACTORY_H_
