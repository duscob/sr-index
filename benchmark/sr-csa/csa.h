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
#include "sr-index/lf.h"
#include "sr-index/construct.h"
#include "sr-index/sequence_ops.h"

#include "index_base.h"

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>>
class CSA : public IndexBaseWithExternalStorage {
 public:

  explicit CSA(std::reference_wrapper<ExternalStorage> t_storage) : IndexBaseWithExternalStorage(t_storage) {}

  void load(sdsl::cache_config t_config) override {
    TSource source(std::ref(t_config));
    loadAllItems(source);
  }

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += serializeItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += serializeItem<TBvMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks");
    written_bytes += serializeRank<TBvMark>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks_rank");
    written_bytes += serializeSelect<TBvMark>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, out, child, "marks_select");

    written_bytes += serializeItem<TMarkToSampleIdx>(
        sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, out, child, "mark_to_sample");

    written_bytes += serializeItem<TSample>(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, out, child, "samples");

    return written_bytes;
  }

 protected:

  virtual void loadAllItems(TSource &t_source) {
    // Create LF function
    auto lf = constructLF(t_source);

    // Create getter for backward search step data
    auto compute_data_backward_search_step = constructComputeDataBackwardSearchStepForPhiForward(t_source);

    // Create function to compute SA values in the range
    auto compute_sa_values = constructComputeSAValues(t_source);

    // Create getter for initial data for backward search step
    auto get_initial_data_backward_search_step = constructGetInitialDataBackwardSearchStep(t_source);

    // Create getter for symbol
    auto get_symbol = constructGetSymbol(t_source);

    auto create_full_range = [](auto tt_seq_size) { return Range{0, tt_seq_size}; };

    auto is_range_empty = [](const auto &tt_range) {
      const auto &[start, end] = tt_range;
      return !(start < end);
    };

    index_.reset(new sri::RIndex(lf,
                                 compute_data_backward_search_step,
                                 compute_sa_values,
                                 n_,
                                 get_initial_data_backward_search_step,
                                 get_symbol,
                                 create_full_range,
                                 is_range_empty));
  }

  struct DataLF {
    std::size_t value = 0;

    struct Run {
      std::size_t start = 0;
      std::size_t end = 0;
      std::size_t rank = 0;
    } run;

    bool operator==(const DataLF &rhs) const {
      return value == rhs.value; // && run.start == rhs.run.start && run.end == rhs.run.end && run.rank == rhs.run.rank;
    }

    bool operator<(const DataLF &rhs) const {
      return value < rhs.value; // || run.start < rhs.run.start || run.end < rhs.run.end || run.rank < rhs.run.rank;
    }
  };

  using Char = sri::uchar;
  using Position = std::size_t;
  struct RangeLF;

  struct Range {
    Position start;
    Position end;

    Range &operator=(const RangeLF &t_range) {
      start = t_range.start.value;
      end = t_range.end.value;
      return *this;
    }
  };

  struct RangeLF {
    DataLF start;
    DataLF end;
  };

  auto constructLF(TSource &t_source) {
    auto cref_alphabet = loadItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, t_source);
    auto cumulative = sri::RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));
    n_ = cumulative[cref_alphabet.get().sigma];

    auto cref_psi_core = loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);
    auto psi_rank = [cref_psi_core](auto tt_c, auto tt_rnk) {
      DataLF data;
      auto report =
          [&data](const auto &tt_rank, const auto &tt_run_start, const auto &tt_run_end, const auto &tt_run_rank) {
            data = DataLF{tt_rank, {tt_run_start, tt_run_end, tt_run_rank}};
          };
      cref_psi_core.get().rank(tt_c, tt_rnk, report);
      return data;
    };

    auto create_range = [](auto tt_c_before_sp, auto tt_c_until_ep, const auto &tt_smaller_c) -> RangeLF {
      tt_c_before_sp.value += tt_smaller_c;
      tt_c_until_ep.value += tt_smaller_c;
      return {tt_c_before_sp, tt_c_until_ep};
    };

    RangeLF empty_range;

    return sri::LF(psi_rank, cumulative, create_range, empty_range);
  }

  auto constructComputeDataBackwardSearchStepForPhiForward(TSource &t_source) {
    auto is_lf_trivial_with_psi = [](const auto &tt_range, const auto &tt_c, const auto &tt_next_range) {
      const auto &[next_start, next_end] = tt_next_range;
      const auto &[start, end] = tt_range;
      return !(next_start < next_end) || (!(start < next_start.run.start) && start < next_start.run.end);
    };

    auto create_data = [](const auto &tt_range, const auto &tt_c, const auto &tt_new_range, const auto &tt_step) {
      const auto &[start, end] = tt_range;
      return DataBackwardSearchStep{tt_c, tt_step, start, end};
    };

    return sri::ComputeDataBackwardSearchStep(is_lf_trivial_with_psi, create_data);
  }

  using Value = std::size_t;
  using TFnReport = std::function<void(std::size_t)>;
  using TFnPhiForRange = std::function<void(const Range &, Value, TFnReport)>;
  virtual TFnPhiForRange constructPhiForRange(TSource &t_source) {
    auto bv_mark_rank = loadBVRank<TBvMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_source, true);
    auto bv_mark_select = loadBVSelect<TBvMark>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_source, true);
    auto successor = sri::CircularSoftSuccessor(bv_mark_rank, bv_mark_select, n_);

    auto cref_mark_to_sample_idx =
        loadItem<TMarkToSampleIdx>(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_source);
    auto get_mark_to_sample_idx = sri::RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, true);

    auto cref_samples = loadItem<TSample>(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto get_sample = sri::RandomAccessForCRefContainer(cref_samples);
    sri::SampleValidatorDefault sample_validator_default;

    auto phi = sri::buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, n_);
    auto phi_for_range = [phi](const auto &t_range, std::size_t t_k, auto t_report) {
      auto k = t_k;
      const auto &[start, end] = t_range;
      for (auto i = start; i < end; ++i) {
        k = phi(k).first;
        t_report(k);
      }
    };

    return phi_for_range;
  }

  using DataBackwardSearchStep = sri::DataBackwardSearchStep<Char>;
  using TFnComputeToehold = std::function<Value(const DataBackwardSearchStep &)>;
  virtual TFnComputeToehold constructComputeToeholdForPhiForward(TSource &t_source) {
    auto cref_psi_core = loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);

    auto cref_samples = loadItem<TSample>(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto get_sa_value_for_bwt_run_start = [cref_psi_core, cref_samples](const auto &tt_pos) {
      auto run = cref_psi_core.get().rankRun(tt_pos);
      return cref_samples.get()[run] + 1;
    };

    return sri::buildComputeToeholdForPhiForward(cref_psi_core, get_sa_value_for_bwt_run_start);
  }

  using TFnComputeSAValues = std::function<void(const Range &, const DataBackwardSearchStep &, TFnReport)>;
  virtual TFnComputeSAValues constructComputeSAValues(TSource &t_source) {
    // Create Phi function using PhiForward
    auto phi_for_range = constructPhiForRange(t_source);

    // Create toehold for phi forward
    auto compute_toehold = constructComputeToeholdForPhiForward(t_source);

    // Create ComputeAllValuesWithPhiForRange
    auto compute_sa_values = sri::buildComputeAllValuesWithPhiForwardForRange(phi_for_range, compute_toehold);

    return compute_sa_values;
  }

  auto constructGetInitialDataBackwardSearchStep(TSource &t_source) {
    auto cref_psi_core = loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);

    return sri::GetInitialDataBackwardSearchStep(cref_psi_core.get().getFirstBWTSymbol(), n_ - 1);
  }

  auto constructGetSymbol(TSource &t_source) {
    auto cref_alphabet = loadItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, t_source);

    auto get_symbol = [cref_alphabet](const auto &tt_c) { return cref_alphabet.get().char2comp[tt_c]; };
    return get_symbol;
  }

};

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBvMark, typename TMarkToSampleIdx, typename TSample>
void construct(CSA<t_width, TAlphabet, TPsiCore, TBvMark, TMarkToSampleIdx, TSample> &t_index,
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
    if (!sdsl::cache_file_exists<TBvMark>(KEY_BWT_RUN_LAST_TEXT_POS, t_config)) {
      sri::constructBitVectorFromIntVector<TBvMark>(KEY_BWT_RUN_LAST_TEXT_POS, t_config, n);
    }
  }

  t_index.load(t_config);
}

#endif //SRI_BENCHMARK_SR_CSA_CSA_H_
