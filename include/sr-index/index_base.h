//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_INDEX_BASE_H_
#define SRI_INDEX_BASE_H_

#include <any>
#include <functional>
#include <map>
#include <string>
#include <variant>

#include <sdsl/io.hpp>

#include "config.h"

namespace sri {

template <typename TSequence = Alphabet<>::string_type>
class CountIndex {
 public:
  virtual ~CountIndex() = default;

  virtual std::pair<std::size_t, std::size_t> Count(const TSequence& _pattern) const = 0;
};

template <typename TSequence = Alphabet<>::string_type>
class LocateIndex {
 public:
  virtual ~LocateIndex() = default;

  virtual std::vector<std::size_t> Locate(const TSequence& _pattern) const = 0;
};

//~~~~~~~


template <typename TSequence,
          typename TBackwardNav,
          typename TGetSymbol,
          typename TCreateFullRange,
          typename TIsRangeEmpty>
class RIndexCountBase : public CountIndex<TSequence> {
 public:
  RIndexCountBase(const TSequence& t_sequence,
                  const TBackwardNav& t_lf,
                  std::size_t t_bwt_size,
                  const TGetSymbol& t_get_symbol,
                  const TCreateFullRange& t_create_full_range,
                  const TIsRangeEmpty& t_is_range_empty)
      : lf_{t_lf},
        bwt_size_{t_bwt_size},
        get_symbol_{t_get_symbol},
        create_full_range_{t_create_full_range},
        is_range_empty_{t_is_range_empty} {}

  std::pair<std::size_t, std::size_t> Count(const TSequence& t_pattern) const override {
    std::pair<std::size_t, std::size_t> range;
    auto report = [&range](const auto& tt_range) {
      const auto& [start, end] = tt_range;
      range = {start, end};
    };

    Count(t_pattern, report);

    return range;
  }

  template <typename TPattern, typename TReport>
  void Count(const TPattern& t_pattern, TReport& t_report) const {
    auto range = create_full_range_(bwt_size_);

    for (auto it = rbegin(t_pattern); it != rend(t_pattern) && !is_range_empty_(range); ++it) {
      auto c = get_symbol_(*it);
      range = lf_(range, c);
    }

    t_report(range);
  }

 protected:
  TBackwardNav lf_;
  std::size_t bwt_size_;
  TGetSymbol get_symbol_;
  TCreateFullRange create_full_range_;
  TIsRangeEmpty is_range_empty_;
};

//~~~~~~~


template <typename TSequence,
          typename TBackwardNav,
          typename TUpdateToeholdData,
          typename TComputeAllValues,
          typename TGetInitialToeholdData,
          typename TGetSymbol,
          typename TCreateFullRange,
          typename TIsRangeEmpty>
class RIndexBase : public RIndexCountBase<TSequence, TBackwardNav, TGetSymbol, TCreateFullRange, TIsRangeEmpty>,
                   public LocateIndex<TSequence> {
 public:
  using CountBase = RIndexCountBase<TSequence, TBackwardNav, TGetSymbol, TCreateFullRange, TIsRangeEmpty>;

  RIndexBase(const TSequence& t_sequence,
             const TBackwardNav& t_lf,
             const TUpdateToeholdData& t_update_toehold_data,
             const TComputeAllValues& t_compute_all_values,
             std::size_t t_bwt_size,
             const TGetInitialToeholdData& t_get_initial_toehold_data,
             const TGetSymbol& t_get_symbol,
             const TCreateFullRange& t_create_full_range,
             const TIsRangeEmpty& t_is_range_empty)
      : CountBase(t_sequence, t_lf, t_bwt_size, t_get_symbol, t_create_full_range, t_is_range_empty),
        update_toehold_data_{t_update_toehold_data},
        compute_all_values_{t_compute_all_values},
        get_initial_toehold_data_{t_get_initial_toehold_data} {}

  std::vector<std::size_t> Locate(const TSequence& t_pattern) const override {
    std::vector<std::size_t> values;
    auto report = [&values](const auto& v) {
      values.emplace_back(v);
    };

    Locate(t_pattern, report);

    return values;
  }

  template <typename TPattern, typename TReport>
  void Locate(const TPattern& t_pattern, TReport& t_report) const {
    auto range = this->create_full_range_(this->bwt_size_);

    auto i = t_pattern.size() - 1;
    // TODO use default value (step == 0) instead of get_initial_toehold_data_
    auto toehold_data = get_initial_toehold_data_(i);

    for (auto it = rbegin(t_pattern); it != rend(t_pattern) && !this->is_range_empty_(range); ++it, --i) {
      auto c = this->get_symbol_(*it);

      auto next_range = this->lf_(range, c);
      update_toehold_data_(range, next_range, c, i, toehold_data);

      range = next_range;
    }

    if (!this->is_range_empty_(range)) {
      compute_all_values_(range, toehold_data, t_report);
    }
  }

 private:
  TUpdateToeholdData update_toehold_data_;
  TComputeAllValues compute_all_values_;
  TGetInitialToeholdData get_initial_toehold_data_;
};

//~~~~~~~


using GenericStorage = std::map<std::string, std::any>;

// Disambiguates an in-memory storage key by the stored type, so a `GenericStorage` shared
// across multiple indexes (e.g. via `ExternalGenericStorage` in the benchmarks) can hold
// different concrete components under the same logical t_key without collision. Used for
// every storage_ access in the helpers below.
template <typename TItem>
inline std::string storageKey(const std::string& t_key) {
  return t_key + "_" + sdsl::util::class_to_hash(TItem{});
}

template <typename TItem>
const TItem* get(const GenericStorage& t_storage, const std::string& t_key) {
  auto it = t_storage.find(t_key);
  return (it != t_storage.end()) ? std::any_cast<TItem>(&it->second) : nullptr;
}

template <typename TItem>
const TItem* set(GenericStorage& t_storage, const std::string& t_key, TItem&& t_item) {
  auto [it, inserted] = t_storage.emplace(t_key, t_item);
  return std::any_cast<TItem>(&it->second);
}

//~~~~~~~


template <typename TStorage = GenericStorage>
class IndexBaseWithExternalStorage {
 public:
  using Storage = TStorage;

  explicit IndexBaseWithExternalStorage(const TStorage& t_storage) : storage_{t_storage} {}

  IndexBaseWithExternalStorage() = default;
  virtual ~IndexBaseWithExternalStorage() = default;

  virtual void load(Config t_config) = 0;
  virtual void load(std::istream& in) = 0;

  typedef std::size_t size_type;

  virtual size_type serialize(std::ostream& out) const {
    return serialize(out, nullptr, "");
  }

  virtual size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const = 0;

 protected:
  enum class ItemKey : unsigned char {
    ALPHABET = 0,
    NAVIGATE,
    MARKS,
    SAMPLES,
    MARK_TO_SAMPLE,
    SAMPLES_IDX,
    RUN_CUMULATIVE_COUNT,
    VALID_MARKS,
    VALID_AREAS,
    NUM_ITEMS
  };

  auto& key(const ItemKey& t_key_enum) {
    return keys_[static_cast<unsigned char>(t_key_enum)];
  }

  const auto& key(const ItemKey& t_key_enum) const {
    return keys_[static_cast<unsigned char>(t_key_enum)];
  }

  using TSource = std::variant<std::reference_wrapper<Config>, std::reference_wrapper<std::istream>>;

  template <typename TItem>
  auto loadRawItem(const std::string& t_key, TSource& t_source, bool t_add_type_hash = false) {
    // The in-memory storage key is always type-discriminated (decoupled from t_add_type_hash,
    // which controls the on-disk filename only). This lets a shared GenericStorage hold
    // distinct components under the same logical t_key for different indexes, and ensures
    // a (t_key, TItem) pair always hits the same cache entry regardless of how subsequent
    // callers spell the t_add_type_hash flag.
    const auto sk = storageKey<TItem>(t_key);
    auto item = get<TItem>(storage_, sk);
    if (!item) {
      // Store a default-constructed item first, then load into it in place. Loading then
      // moving into storage would invalidate SDSL rank/select support pointers (e.g.
      // rank_support_sd::m_v) wired up via set_vector during deserialization.
      auto* mutable_item = const_cast<TItem*>(set(storage_, sk, TItem{}));
      load(*mutable_item, t_source, t_key, t_add_type_hash);
      item = mutable_item;
    }
    return item;
  }

  template <typename TItem>
  auto load(TItem& t_item, TSource t_source, const std::string& t_key, bool t_add_type_hash) {
    std::visit(
        [&t_item, &t_key, &t_add_type_hash, this](auto&& tt_source) {
          return this->load(t_item, tt_source.get(), t_key, t_add_type_hash);
        },
        t_source);
  }

  template <typename TItem>
  auto load(TItem& t_item, const sdsl::cache_config& t_config, const std::string& t_key, bool t_add_type_hash) {
    if (!sdsl::load_from_cache(t_item, t_key, t_config, t_add_type_hash))
      throw std::invalid_argument("File not found (Key: '" + t_key + "')");
  }

  template <typename TItem>
  auto load(TItem& t_item, std::istream& t_in, const std::string& t_key, bool) {
    sdsl::load(t_item, t_in);
  }

  template <typename TItem>
  auto loadItem(const std::string& t_key, TSource& t_source, bool t_add_type_hash = false) {
    auto item = loadRawItem<TItem>(t_key, t_source, t_add_type_hash);
    return std::cref(*item);
  }

  template <typename TBv, typename TBvRank = typename TBv::rank_1_type>
  auto loadBVRank(const std::string& t_key, TSource& t_source, bool t_add_type_hash = false) {
    // The rank's storage entry is type-discriminated by TBvRank — different TBv variants
    // sharing t_key produce different rank-support types, and must not collide on the
    // shared storage_ map.
    const auto sk_rank = storageKey<TBvRank>(t_key + "_rank");
    auto item_rank = get<TBvRank>(storage_, sk_rank);
    if (!item_rank) {
      auto item_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      auto* mutable_rank = const_cast<TBvRank*>(set(storage_, sk_rank, TBvRank{}));
      load(*mutable_rank, t_source, t_key, t_add_type_hash);
      mutable_rank->set_vector(item_bv);
      item_rank = mutable_rank;
    }

    return std::cref(*item_rank);
  }

  template <typename TBv, typename TBvSelect = typename TBv::select_1_type>
  auto loadBVSelect(const std::string& t_key, TSource& t_source, bool t_add_type_hash = false) {
    const auto sk_select = storageKey<TBvSelect>(t_key + "_select");
    auto item_select = get<TBvSelect>(storage_, sk_select);
    if (!item_select) {
      auto item_bv = loadRawItem<TBv>(t_key, t_source, t_add_type_hash);

      auto* mutable_select = const_cast<TBvSelect*>(set(storage_, sk_select, TBvSelect{}));
      load(*mutable_select, t_source, t_key, t_add_type_hash);
      mutable_select->set_vector(item_bv);
      item_select = mutable_select;
    }

    return std::cref(*item_select);
  }

  template <typename TItem>
  std::size_t serializeItem(const std::string& t_key,
                            std::ostream& out,
                            sdsl::structure_tree_node* v,
                            const std::string& name) const {
    // Mirror the type-discriminated storage key used by loadRawItem.
    auto item = get<TItem>(storage_, storageKey<TItem>(t_key));
    if (item) {
      return sdsl::serialize(*item, out, v, name);
    }
    return sdsl::serialize_empty_object<TItem>(out, v, name);
  }

  template <typename TItem, typename TItemRank = typename TItem::rank_1_type>
  std::size_t serializeRank(const std::string& t_key,
                            std::ostream& out,
                            sdsl::structure_tree_node* v,
                            const std::string& name) const {
    return serializeItem<TItemRank>(t_key + "_rank", out, v, name);
  }

  template <typename TItem, typename TItemSelect = typename TItem::select_1_type>
  std::size_t serializeSelect(const std::string& t_key,
                              std::ostream& out,
                              sdsl::structure_tree_node* v,
                              const std::string& name) const {
    return serializeItem<TItemSelect>(t_key + "_select", out, v, name);
  }

  //********************
  //********************
  //********************

  std::size_t n_ = 0;
  TStorage storage_;
  std::array<std::string, static_cast<u_int8_t>(ItemKey::NUM_ITEMS)> keys_;
};

//~~~~~~~


template <typename TStorage = GenericStorage, typename TSequence = Alphabet<>::string_type>
class CountIndexExtStorage : public CountIndex<TSequence>, public IndexBaseWithExternalStorage<TStorage> {
 public:
  explicit CountIndexExtStorage(const TStorage& t_storage) : IndexBaseWithExternalStorage<TStorage>(t_storage) {}

  CountIndexExtStorage() = default;

  std::pair<std::size_t, std::size_t> Count(const TSequence& t_pattern) const override {
    return index_->Count(t_pattern);
  }

 protected:
  std::shared_ptr<CountIndex<TSequence>> index_ = nullptr;
};

//~~~~~~~


template <typename TStorage = GenericStorage, typename TSequence = Alphabet<>::string_type>
class LocateIndexExtStorage : public LocateIndex<TSequence>, public IndexBaseWithExternalStorage<TStorage> {
 public:
  explicit LocateIndexExtStorage(const TStorage& t_storage) : IndexBaseWithExternalStorage<TStorage>(t_storage) {}

  LocateIndexExtStorage() = default;

  std::vector<std::size_t> Locate(const TSequence& t_pattern) const override {
    return index_->Locate(t_pattern);
  }

 protected:
  std::shared_ptr<LocateIndex<TSequence>> index_ = nullptr;
};

//~~~~~~~


template <typename TSequence,
          typename TBackwardNav,
          typename TGetLastValue,
          typename TComputeAllValues,
          typename TGetFinalValue,
          typename TGetSymbol>
auto buildSharedPtrRIndex(const TBackwardNav& t_lf,
                          const TGetLastValue& t_get_last_value,
                          const TComputeAllValues& t_compute_all_values,
                          std::size_t t_bwt_size,
                          const TGetFinalValue& t_get_final_sa_value,
                          const TGetSymbol& t_get_symbol) {
  using Range = std::pair<std::size_t, std::size_t>;
  using TFnCreateFullRange = std::function<Range(std::size_t)>;
  auto create_full_range = [](auto tt_seq_size) {
    return Range{0, tt_seq_size - 1};
  };

  using TFnIsRangeEmpty = std::function<bool(const Range&)>;
  auto is_range_empty = [](const auto& tt_range) {
    return tt_range.second < tt_range.first;
  };

  return std::make_shared<RIndexBase<TSequence,
                                     TBackwardNav,
                                     TGetLastValue,
                                     TComputeAllValues,
                                     TGetFinalValue,
                                     TGetSymbol,
                                     TFnCreateFullRange,
                                     TFnIsRangeEmpty>>(TSequence{},
                                                       t_lf,
                                                       t_get_last_value,
                                                       t_compute_all_values,
                                                       t_bwt_size,
                                                       t_get_final_sa_value,
                                                       t_get_symbol,
                                                       create_full_range,
                                                       is_range_empty);
}

}  // namespace sri

#endif  // SRI_INDEX_BASE_H_
