//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/25/21.
//

#ifndef SRI_BENCHMARK_SR_INDEX_FACTORY_H_
#define SRI_BENCHMARK_SR_INDEX_FACTORY_H_

#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"
#include "config.h"

using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

template<uint8_t t_width = 8>
class Factory {
 public:
  enum class IndexEnum {
    R_INDEX = 0,
    SR_INDEX,
    SR_INDEX_VM,
    SR_INDEX_VA,
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
      case IndexEnum::R_INDEX: {
        auto idx = std::make_shared<sri::RIndex<ExternalGenericStorage>>(std::ref(storage_));
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_INDEX: {
        auto idx = std::make_shared<sri::SrIndex<ExternalGenericStorage>>(std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_INDEX_VM: {
        auto idx = std::make_shared<sri::SrIndexValidMark<ExternalGenericStorage>>(
            std::ref(storage_), t_config.sampling_size);
        idx->load(config_);
        index = {idx, sdsl::size_in_bytes(*idx)};
        break;
      }

      case IndexEnum::SR_INDEX_VA: {
        auto idx = std::make_shared<sri::SrIndexValidArea<ExternalGenericStorage>>(
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

#endif //SRI_BENCHMARK_SR_INDEX_FACTORY_H_
