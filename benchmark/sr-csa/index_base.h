//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_INDEX_BASE_H_
#define SRI_BENCHMARK_SR_CSA_INDEX_BASE_H_

#include <map>
#include <string>
#include <any>
#include <functional>
#include <variant>

#include <sdsl/io.hpp>

#include "sr-index/r_index.h"

using ExternalStorage = std::map<std::string, std::any>;

enum class SrIndexKey : unsigned char {
  ALPHABET = 0,
  NAVIGATE,
  MARKS,
  SAMPLES,
  MARK_TO_SAMPLE,
  SAMPLES_IDX,
  RUN_CUMULATIVE_COUNT,
  VALID_MARKS,
  VALID_AREAS
};

auto toInt(SrIndexKey t_k) {
  return static_cast<unsigned char>(t_k);
}

class IndexBaseWithExternalStorage : public sri::LocateIndex {
 public:
  explicit IndexBaseWithExternalStorage(const std::reference_wrapper<ExternalStorage> &t_storage)
      : storage_{t_storage} {
  }

  std::vector<std::size_t> Locate(const std::string &t_pattern) const override {
    return index_->Locate(t_pattern);
  }

  auto sizeSequence() const { return n_; }

  virtual void load(sdsl::cache_config t_config) = 0;

  typedef std::size_t size_type;
//  virtual void load(std::istream &in) = 0;

  virtual size_type serialize(std::ostream &out) const {
    return serialize(out, nullptr, "");
  }

  virtual size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const = 0;

 protected:

  auto &key(const SrIndexKey &t_key_enum) {
    return keys_[toInt(t_key_enum)];
  }

  const auto &key(const SrIndexKey &t_key_enum) const {
    return keys_[toInt(t_key_enum)];
  }

  using TSource = std::variant<std::reference_wrapper<sdsl::cache_config>, std::reference_wrapper<std::istream>>;

  template<typename TItem>
  auto loadRawItem(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto it = storage_.get().find(t_key);
    if (it == storage_.get().end()) {
      TItem data;
      load(data, t_source, t_key, t_add_type_hash);

      bool inserted = false;
      std::tie(it, inserted) = storage_.get().emplace(t_key, std::move(data));
    }

    return it;
  }

  template<typename TItem>
  auto load(TItem &t_item, TSource t_source, const std::string &t_key, bool t_add_type_hash) {
    std::visit([&t_item, &t_key, &t_add_type_hash, this](auto &&tt_source) {
                 return load(t_item, tt_source.get(), t_key, t_add_type_hash);
               },
               t_source);
  }

  template<typename TItem>
  auto load(TItem &t_item, const sdsl::cache_config &t_config, const std::string &t_key, bool t_add_type_hash) {
    sdsl::load_from_cache(t_item, t_key, t_config, t_add_type_hash);
  }

  template<typename TItem>
  auto load(TItem &t_item, std::istream &t_in, const std::string &t_key, bool) {
    sdsl::load(t_item, t_in);
  }

  template<typename TItem>
  auto loadItem(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto it = loadRawItem<TItem>(t_key, t_source, t_add_type_hash);
    return std::cref(std::any_cast<const TItem &>(it->second));
  }

  template<typename TBv, typename TBvRank = typename TBv::rank_1_type>
  auto loadBVRank(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_rank = t_key + "_rank";
    auto it = storage_.get().find(key_rank);
    if (it == storage_.get().end()) {
      auto it_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      TBvRank rank;
      load(rank, t_source, t_key, t_add_type_hash);
      rank.set_vector(std::any_cast<TBv>(&it_bv->second));

      bool inserted = false;
      std::tie(it, inserted) = storage_.get().emplace(key_rank, std::move(rank));
    }

    return std::cref(std::any_cast<const TBvRank &>(it->second));
  }

  template<typename TBv, typename TBvSelect = typename TBv::select_1_type>
  auto loadBVSelect(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_select = t_key + "_select";
    auto it = storage_.get().find(key_select);
    if (it == storage_.get().end()) {
      auto it_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      TBvSelect select;
      load(select, t_source, t_key, t_add_type_hash);
      select.set_vector(std::any_cast<TBv>(&it_bv->second));

      bool inserted = false;
      std::tie(it, inserted) = storage_.get().emplace(key_select, std::move(select));
    }

    return std::cref(std::any_cast<const TBvSelect &>(it->second));
  }

  template<typename TItem>
  std::size_t serializeItem(const std::string &t_key,
                            std::ostream &out,
                            sdsl::structure_tree_node *v,
                            const std::string &name) const {
    auto it = storage_.get().find(t_key);
    if (it != storage_.get().end()) {
      auto cref_item = std::cref(std::any_cast<const TItem &>(it->second));
      return sdsl::serialize(cref_item.get(), out, v, name);
    }
    return sdsl::serialize_empty_object<TItem>(out, v, name);
  }

  template<typename TItem, typename TItemRank = typename TItem::rank_1_type>
  std::size_t serializeRank(const std::string &t_key,
                            std::ostream &out,
                            sdsl::structure_tree_node *v,
                            const std::string &name) const {
    return serializeItem<TItemRank>(t_key + "_rank", out, v, name);
  }

  template<typename TItem, typename TItemSelect = typename TItem::select_1_type>
  std::size_t serializeSelect(const std::string &t_key,
                              std::ostream &out,
                              sdsl::structure_tree_node *v,
                              const std::string &name) const {
    return serializeItem<TItemSelect>(t_key + "_select", out, v, name);
  }

  //********************
  //********************
  //********************

  std::size_t n_ = 0;
  std::reference_wrapper<ExternalStorage> storage_;
  std::vector<std::string> keys_;

  std::shared_ptr<sri::LocateIndex> index_ = nullptr;
};

#endif //SRI_BENCHMARK_SR_CSA_INDEX_BASE_H_
