//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_FACTORY_H_
#define SRI_BENCHMARK_SR_CSA_FACTORY_H_

#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "csa.h"
#include "sr_csa.h"
#include "csa_raw.h"

template<uint8_t t_width = 8>
class Factory {
 public:
  enum class IndexEnum {
    R_CSA = 0,
    CSA_RAW = 1,
    SR_CSA = 2,
    SR_CSA_SLIM = 3,
    SR_CSA_BV = 4
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
        auto idx = std::make_shared<CSA<t_width>>(std::ref(storage_));
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::CSA_RAW: {
        auto idx = std::make_shared<CSARaw<t_width>>(std::ref(storage_));
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA: {
        auto idx = std::make_shared<SrCSA<t_width>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_SLIM: {
        auto idx = std::make_shared<SrCSASlim<t_width>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }

      case IndexEnum::SR_CSA_BV: {
        auto idx = std::make_shared<SrCSAWithBv<t_width>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        return {idx, sdsl::size_in_bytes(*idx)};
      }
    }

    return {nullptr, 0};
  }

 private:
  sdsl::cache_config config_;

  std::size_t n_ = 0;

  ExternalStorage storage_;
};

#endif //SRI_BENCHMARK_SR_CSA_FACTORY_H_
