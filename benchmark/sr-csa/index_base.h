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

template<typename TItem>
const TItem *get(const std::reference_wrapper<ExternalStorage> &t_storage, const std::string &t_key) {
  auto it = t_storage.get().find(t_key);
  return (it != t_storage.get().end()) ? std::any_cast<TItem>(&it->second) : nullptr;
}

template<typename TItem>
const TItem *set(std::reference_wrapper<ExternalStorage> t_storage, const std::string &t_key, TItem &&t_item) {
  auto[it, inserted] = t_storage.get().emplace(t_key, t_item);
  return std::any_cast<TItem>(&it->second);
}

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

  IndexBaseWithExternalStorage() = default;

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
    auto item = get<TItem>(storage_, t_key);
    if (!item) {
      TItem data;
      load(data, t_source, t_key, t_add_type_hash);

      item = set(storage_, t_key, std::move(data));
    }
    return item;
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
    auto item = loadRawItem<TItem>(t_key, t_source, t_add_type_hash);
    return std::cref(*item);
  }

  template<typename TBv, typename TBvRank = typename TBv::rank_1_type>
  auto loadBVRank(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_rank = t_key + "_rank";
    auto item_rank = get<TBvRank>(storage_, key_rank);
    if (!item_rank) {
      auto item_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      TBvRank rank;
      load(rank, t_source, t_key, t_add_type_hash);
      rank.set_vector(item_bv);

      item_rank = set(storage_, key_rank, std::move(rank));
    }

    return std::cref(*item_rank);
  }

  template<typename TBv, typename TBvSelect = typename TBv::select_1_type>
  auto loadBVSelect(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_select = t_key + "_select";
    auto item_select = get<TBvSelect>(storage_, key_select);
    if (!item_select) {
      auto item_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      TBvSelect select;
      load(select, t_source, t_key, t_add_type_hash);
      select.set_vector(item_bv);

      item_select = set(storage_, key_select, std::move(select));
    }

    return std::cref(*item_select);
  }

  template<typename TItem>
  std::size_t serializeItem(const std::string &t_key,
                            std::ostream &out,
                            sdsl::structure_tree_node *v,
                            const std::string &name) const {
    auto item = get<TItem>(storage_, t_key);
    if (item) {
      return sdsl::serialize(*item, out, v, name);
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
