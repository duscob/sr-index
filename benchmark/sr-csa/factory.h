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
    SR_CSA_SLIM_VA
  };

  struct Config {
    IndexEnum index_type;
    std::size_t sampling_size;
  };

  explicit Factory(sdsl::cache_config t_config) : config_{std::move(t_config)} {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, config_));
    n_ = buf.size();
  }

  auto sizeSequence() const { return n_; }

  std::pair<std::shared_ptr<sri::LocateIndex>, std::size_t> make(const Config &t_config) {
    switch (t_config.index_type) {
      case IndexEnum::R_CSA: {
        auto idx = std::make_shared<sri::RCSA<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::CSA_RAW: {
        auto idx = std::make_shared<sri::CSARaw<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA: {
        auto idx = std::make_shared<sri::SrCSA<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_VM: {
        auto idx = std::make_shared<sri::SrCSAValidMark<sri::SrCSA<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_VA: {
        auto idx = std::make_shared<sri::SrCSAValidArea<sri::SrCSA<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_SLIM: {
        auto idx = std::make_shared<sri::SrCSASlim<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_SLIM_VM: {
        auto idx = std::make_shared<sri::SrCSAValidMark<sri::SrCSASlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_SLIM_VA: {
        auto idx = std::make_shared<sri::SrCSAValidArea<sri::SrCSASlim<ExternalGenericStorage>>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }
    }

    return {nullptr, 0};
  }

 private:
  sdsl::cache_config config_;

  std::size_t n_ = 0;

  sri::GenericStorage storage_;
};

#endif //SRI_BENCHMARK_SR_CSA_FACTORY_H_
