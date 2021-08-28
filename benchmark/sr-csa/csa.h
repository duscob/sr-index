//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/15/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_CSA_H_
#define SRI_BENCHMARK_SR_CSA_CSA_H_

#include <map>
#include <string>
#include <any>
#include <functional>
#include <memory>

#include <sdsl/csa_alphabet_strategy.hpp>

#include "sr-index/r_index.h"
#include "sr-index/psi.h"
#include "sr-index/construct.h"
#include "sr-index/sequence_ops.h"

using ExternalStorage = std::map<std::string, std::any>;

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBVMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>>
class CSA : public sri::LocateIndex {
 public:
  explicit CSA(std::reference_wrapper<ExternalStorage> t_storage) : storage_{t_storage} {}

  std::vector<std::size_t> Locate(const std::string &t_pattern) const override {
    return index_->Locate(t_pattern);
  }

  auto sizeSequence() const { return n_; }

  typedef std::size_t size_type;

  void load(sdsl::cache_config t_config) {
    loadAllItems(t_config);
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, const std::string &name = "") const {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += serializeItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += serializeItem<TBVMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks");
    written_bytes += serializeRank<TBVMark>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks_rank");
    written_bytes += serializeSelect<TBVMark>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks_select");

    written_bytes += serializeItem<TMarkToSampleIdx>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, out, child, "mark_to_sample");

    written_bytes += serializeItem<TSample>(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, out, child, "samples");

    return written_bytes;
  }

 private:

  template<typename TSource>
  void loadAllItems(TSource &t_source) {
    // Create LF function
    auto cref_alphabet = loadItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, t_source);
    auto cumulative = sri::RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));
    n_ = cumulative[cref_alphabet.get().sigma];

    auto cref_psi_core = loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);
    auto psi_rank = [cref_psi_core](auto tt_c, auto tt_rnk) { return cref_psi_core.get().rank(tt_c, tt_rnk); };

    auto lf = sri::LFOnPsi(psi_rank, cumulative);

    // Create getter for backward search step data
    auto compute_data_backward_search_step = sri::buildComputeDataBackwardSearchStepForPhiForward(cref_psi_core);

    // Create Phi function using PhiForward
    auto bv_mark_rank = loadBVRank<TBVMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_source, true);
    auto bv_mark_select = loadBVSelect<TBVMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_source, true);
    auto successor = sri::CircularSoftSuccessor(bv_mark_rank, bv_mark_select, n_);

    auto cref_mark_to_sample_idx =
        loadItem<TMarkToSampleIdx>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_source);
    auto get_mark_to_sample_idx = sri::RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, true);

    auto cref_samples = loadItem<TSample>(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto get_sample = sri::RandomAccessForCRefContainer(cref_samples);
    sri::SampleValidatorDefault sample_validator_default;

    auto phi = sri::buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, n_);
    auto phi_for_range = [phi](const auto &t_range, std::experimental::optional<std::size_t> t_k, auto &t_report) {
      auto k = *t_k;
      for (auto i = t_range.first; i <= t_range.second; ++i) {
        k = phi(k).first;
        t_report(k);
      }
    };

    // Create sample value from position in SA
    auto get_sa_value_for_bwt_run_start = [cref_psi_core, cref_samples](const auto &tt_pos) {
      auto run = cref_psi_core.get().rankRun(tt_pos);
      return cref_samples.get()[run] + 1;
    };

    auto compute_toehold = sri::buildComputeToeholdForPhiForward(cref_psi_core, get_sa_value_for_bwt_run_start);

    // Create ComputeAllValuesWithPhiForRange
    auto compute_all_values = sri::buildComputeAllValuesWithPhiForwardForRange(phi_for_range, compute_toehold);

    // Create getter for initial data for backward search step
    auto get_initial_data_backward_search_step = sri::GetInitialDataBackwardSearchStep(
        cref_psi_core.get().getFirstBWTSymbol(), n_ - 1);

    // Create getter for symbol
    auto get_symbol = [cref_alphabet](const auto &tt_c) { return cref_alphabet.get().char2comp[tt_c]; };
    index_.reset(new sri::RIndex(lf,
                                 compute_data_backward_search_step,
                                 compute_all_values,
                                 n_,
                                 get_initial_data_backward_search_step,
                                 get_symbol));
  }

  template<typename TItem, typename TSource>
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
  auto load(TItem &t_item, const sdsl::cache_config &t_config, const std::string &t_key, bool t_add_type_hash) {
    sdsl::load_from_cache(t_item, t_key, t_config, t_add_type_hash);
  }

  template<typename TItem>
  auto load(TItem &t_item, std::istream &t_in, const std::string &t_key, bool) {
    sdsl::load(t_item, t_in);
  }

  template<typename TItem, typename TSource>
  auto loadItem(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto it = loadRawItem<TItem>(t_key, t_source, t_add_type_hash);
    return std::cref(std::any_cast<const TItem &>(it->second));
  }

  template<typename TBV, typename TSource>
  auto loadBVRank(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_rank = t_key + "_rank";
    auto it = storage_.get().find(key_rank);
    if (it == storage_.get().end()) {
      auto it_bv = loadRawItem<TBV>(t_key, t_source, t_add_type_hash);

      typename TBV::rank_1_type rank;
      load(rank, t_source, t_key, t_add_type_hash);
      rank.set_vector(std::any_cast<TBV>(&it_bv->second));

      bool inserted = false;
      std::tie(it, inserted) = storage_.get().emplace(key_rank, std::move(rank));
    }

    return std::cref(std::any_cast<const typename TBV::rank_1_type &>(it->second));
  }

  template<typename TBV, typename TSource>
  auto loadBVSelect(const std::string &t_key, TSource &t_source, bool t_add_type_hash = false) {
    auto key_select = t_key + "_select";
    auto it = storage_.get().find(key_select);
    if (it == storage_.get().end()) {
      auto it_bv = loadRawItem<TBV>(t_key, t_source, t_add_type_hash);

      typename TBV::select_1_type select;
      load(select, t_source, t_key, t_add_type_hash);
      select.set_vector(std::any_cast<TBV>(&it_bv->second));

      bool inserted = false;
      std::tie(it, inserted) = storage_.get().emplace(key_select, std::move(select));
    }

    return std::cref(std::any_cast<const typename TBV::select_1_type &>(it->second));
  }

  template<typename TItem>
  size_type serializeItem(const std::string &t_key,
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

  template<typename TItem>
  size_type serializeRank(const std::string &t_key,
                          std::ostream &out,
                          sdsl::structure_tree_node *v,
                          const std::string &name) const {
    return serializeItem<typename TItem::rank_1_type>(t_key + "_rank", out, v, name);
  }

  template<typename TItem>
  size_type serializeSelect(const std::string &t_key,
                            std::ostream &out,
                            sdsl::structure_tree_node *v,
                            const std::string &name) const {
    return serializeItem<typename TItem::select_1_type>(t_key + "_select", out, v, name);
  }

  //********************
  //********************
  //********************

  std::size_t n_ = 0;
  std::reference_wrapper<ExternalStorage> storage_;

  std::unique_ptr<sri::LocateIndex> index_;
};

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBVMark, typename TMarkToSampleIdx, typename TSample>
void construct(CSA<t_width, TAlphabet, TPsiCore, TBVMark, TMarkToSampleIdx, TSample> &t_index,
               sdsl::cache_config &t_config) {
  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = bwt_buf.size();
  }

  {
    // Construct Successor on the text positions of BWT run last letter
    auto event = sdsl::memory_monitor::event("Successor");
    const auto KEY_BWT_RUN_LAST_TEXT_POS = sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS;
    if (!sdsl::cache_file_exists<TBVMark>(KEY_BWT_RUN_LAST_TEXT_POS, t_config)) {
      sri::constructBitVectorFromIntVector<TBVMark>(KEY_BWT_RUN_LAST_TEXT_POS, t_config, n);
    }
  }

  t_index.load(t_config);
}

#endif //SRI_BENCHMARK_SR_CSA_CSA_H_
