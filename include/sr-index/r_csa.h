//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/15/21.
//

#ifndef SRI_R_CSA_H_
#define SRI_R_CSA_H_

#include <any>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include <sdsl/csa_alphabet_strategy.hpp>

#include "alphabet.h"
#include "config.h"
#include "construct.h"
#include "index_base.h"
#include "lf.h"
#include "psi.h"
#include "sequence_ops.h"

namespace sri {

template <typename TStorage = GenericStorage,
          typename TAlphabet = Alphabet<>,
          typename TPsiRLE = PsiCoreRLE<>,
          typename TBvMark = sdsl::sd_vector<>,
          typename TMarkToSampleIdx = sdsl::int_vector<>,
          typename TSample = sdsl::int_vector<>>
class RCSA : public LocateIndexExtStorage<TStorage, typename TAlphabet::string_type> {
 public:
  using Alphabet = TAlphabet;
  using Samples = TSample;
  using BvMarks = TBvMark;
  using MarksToSamples = TMarkToSampleIdx;
  using Base = LocateIndexExtStorage<TStorage, typename TAlphabet::string_type>;

  explicit RCSA(const TStorage& t_storage) : Base(t_storage) {}

  RCSA() = default;

  virtual ~RCSA() = default;

  void load(Config t_config) override {
    TSource source(std::ref(t_config));
    loadInner(source, t_config.keys);
  }

  void load(std::istream& in) override {
    TSource source(std::ref(in));
    loadInner(source, createDefaultKeys<TAlphabet::int_width>());
  }

  using typename Base::ItemKey;
  using typename Base::size_type;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(key(ItemKey::ALPHABET), out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(key(ItemKey::NAVIGATE), out, child, "psi");

    written_bytes += this->template serializeItem<TSample>(key(ItemKey::SAMPLES), out, child, "samples");

    written_bytes += this->template serializeItem<TBvMark>(key(ItemKey::MARKS), out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(key(ItemKey::MARKS), out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(key(ItemKey::MARKS), out, child, "marks_select");

    written_bytes +=
        this->template serializeItem<TMarkToSampleIdx>(key(ItemKey::MARK_TO_SAMPLE), out, child, "mark_to_sample");

    return written_bytes;
  }

 protected:
  using typename Base::TSource;

  virtual void loadInner(TSource& t_source, const JSON& t_keys) {
    setupKeyNames(t_keys);
    loadAllItems(t_source);
    constructIndex(t_source);
  }

  using Base::key;

  virtual void setupKeyNames(const JSON& t_keys) {
    using namespace sri::conf;
    key(ItemKey::ALPHABET) = t_keys[kAlphabet];
    key(ItemKey::NAVIGATE) = t_keys[kPsi][kBase];
    key(ItemKey::SAMPLES) = t_keys[kPsi][kHead][kTextPos];
    key(ItemKey::MARKS) = t_keys[kPsi][kTail][kTextPos];
    key(ItemKey::MARK_TO_SAMPLE) = t_keys[kPsi][kTail][kTextPosAsc][kLink];
  }

  virtual void loadAllItems(TSource& t_source) {
    this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);

    this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);

    this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source, true);

    this->template loadItem<TBvMark>(key(ItemKey::MARKS), t_source, true);
    this->template loadBVRank<TBvMark>(key(ItemKey::MARKS), t_source, true);
    this->template loadBVSelect<TBvMark>(key(ItemKey::MARKS), t_source, true);

    this->template loadItem<TMarkToSampleIdx>(key(ItemKey::MARK_TO_SAMPLE), t_source, true);
  }

  virtual void constructIndex(TSource& t_source) {
    this->index_.reset(
        new RIndexBase{typename TAlphabet::string_type{},  //
                       constructLF(t_source),              //
                       constructComputeDataBackwardSearchStep(
                           [](const Range& tt_range, auto tt_c, const RangeLF& tt_next_range, std::size_t tt_step) {
                             const auto& [start, end] = tt_next_range;
                             return DataBackwardSearchStep{tt_step, RunData{tt_c, start.run.rank}};
                           }),                                                                                       //
                       constructComputeSAValues(constructPhiForRange(t_source), constructComputeToehold(t_source)),  //
                       this->n_,                                                                                     //
                       [](const auto& tt_step) {
                         return DataBackwardSearchStep{0, RunData{0, 0}};
                       },                             //
                       constructGetSymbol(t_source),  //
                       [](auto tt_seq_size) {
                         return Range{0, tt_seq_size};
                       },  //
                       constructIsRangeEmpty()});
  }

  struct DataLF {
    std::size_t value = 0;

    struct Run {
      std::size_t start = 0;
      std::size_t end = 0;
      std::size_t rank = 0;
    } run;

    bool operator==(const DataLF& rhs) const {
      return value == rhs.value;
      // && run.start == rhs.run.start && run.end == rhs.run.end && run.rank == rhs.run.rank;
    }

    bool operator<(const DataLF& rhs) const {
      return value < rhs.value;
      // || run.start < rhs.run.start || run.end < rhs.run.end || run.rank < rhs.run.rank;
    }
  };

  using Char = typename TAlphabet::comp_char_type;
  using Position = std::size_t;
  struct RangeLF;

  struct Range {
    Position start;
    Position end;

    Range& operator=(const RangeLF& t_range) {
      start = t_range.start.value;
      end = t_range.end.value;
      return *this;
    }
  };

  struct RangeLF {
    DataLF start;
    DataLF end;
  };

  auto constructIsRangeEmpty() {
    return [](const Range& tt_range) {
      return !(tt_range.start < tt_range.end);
    };
  }

  auto constructLF(TSource& t_source) {
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);
    auto cumulative = RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));
    this->n_ = cumulative[cref_alphabet.get().sigma];

    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);
    auto psi_rank = [cref_psi_core](auto tt_c, auto tt_rnk) {
      DataLF data;
      auto report =
          [&data](const auto& tt_rank, const auto& tt_run_start, const auto& tt_run_end, const auto& tt_run_rank) {
            data = DataLF{tt_rank, {tt_run_start, tt_run_end, tt_run_rank}};
          };
      cref_psi_core.get().rank(tt_c, tt_rnk, report);
      return data;
    };

    auto create_range = [](auto tt_c_before_sp, auto tt_c_until_ep, const auto& tt_smaller_c) -> RangeLF {
      tt_c_before_sp.value += tt_smaller_c;
      tt_c_until_ep.value += tt_smaller_c;
      return {tt_c_before_sp, tt_c_until_ep};
    };

    RangeLF empty_range;

    return LF(psi_rank, cumulative, create_range, empty_range);
  }

  struct RunData {
    Char c;                    // Character for LF step in the range
    std::size_t partial_rank;  // Rank of first run item (total rank for symbol c and partial rank for all symbols)
  };

  using DataBackwardSearchStep = sri::DataBackwardSearchStep<RunData>;

  template <typename TCreateData>
  auto constructComputeDataBackwardSearchStep(const TCreateData& t_create_data) {
    return buildComputeDataBackwardSearchStepForPhiForward(t_create_data);
  }

  using Value = std::size_t;

  auto constructPhiForRange(TSource& t_source) {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(key(ItemKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(key(ItemKey::MARKS), t_source, true);
    auto successor = CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(key(ItemKey::MARK_TO_SAMPLE), t_source);
    auto get_mark_to_sample_idx = RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, true);

    auto cref_samples = this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source);
    auto get_sample = RandomAccessForCRefContainer(cref_samples);
    SampleValidatorDefault sample_validator_default;

    auto phi = buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);
    auto phi_for_range = [phi](const auto& t_range, std::size_t t_k, auto t_report) {
      auto k = t_k;
      const auto& [start, end] = t_range;
      for (auto i = start; i < end; ++i) {
        k = phi(k).first;
        t_report(k);
      }
    };

    return phi_for_range;
  }

  auto constructComputeToehold(TSource& t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);

    auto cref_samples = this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source);
    auto get_sa_value_for_bwt_run_start = [cref_psi_core, cref_samples](const RunData& tt_run_data) {
      auto n_prev_runs = cref_psi_core.get().rankCharRun(tt_run_data.c);
      return cref_samples.get()[n_prev_runs + tt_run_data.partial_rank] + 1;
    };

    return buildComputeToeholdForPhiForward(get_sa_value_for_bwt_run_start, cref_psi_core.get().size());
  }

  template <typename TPhiRange, typename TComputeToehold>
  auto constructComputeSAValues(const TPhiRange& t_phi_range, const TComputeToehold& t_compute_toehold) {
    auto update_range = [](Range tt_range) {
      auto& [start, end] = tt_range;
      ++start;
      return tt_range;
    };

    return ComputeAllValuesWithPhiForRange(t_phi_range, t_compute_toehold, update_range);
  }

  auto constructGetSymbol(TSource& t_source) {
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);

    auto get_symbol = [cref_alphabet](auto tt_c) {
      return cref_alphabet.get().char2comp[tt_c];
    };
    return get_symbol;
  }
};

template <typename... TArgs>
void constructItems(RCSA<TArgs...>& t_index, Config& t_config) {
  using Index = RCSA<TArgs...>;
  using namespace sri::conf;
  constexpr auto width = Index::Alphabet::int_width;
  const auto& keys = t_config.keys;

  constructIndexBaseItems<width>(t_config.data_path, t_config);

  // Construct Psi
  if (!cache_file_exists(keys[kPsi][kBase], t_config)) {
    auto event = sdsl::memory_monitor::event("Psi");
    constructPsi<width>(t_config);
  }

  // Construct Psi Runs
  if (!sdsl::cache_file_exists<sdsl::int_vector<>>(keys[kPsi][kHead][kTextPos], t_config)) {
    auto event = sdsl::memory_monitor::event("Psi Runs");
    constructPsiRuns<width>(t_config);
  }

  // Construct Samples for the template type
  if (!std::is_same_v<typename Index::Samples, sdsl::int_vector<>>) {
    auto event = sdsl::memory_monitor::event("Samples");
    sdsl::int_vector<> samples_iv;
    sdsl::load_from_cache(samples_iv, keys[kPsi][kHead][kTextPos], t_config, true);

    auto samples = sri::construct<typename Index::Samples>(samples_iv);
    sri::store_to_cache(samples, keys[kPsi][kHead][kTextPos], t_config, true);
  }

  // Construct Successor on the text positions of Psi run last item
  if (!sdsl::cache_file_exists<typename Index::BvMarks>(keys[kPsi][kTail][kTextPos], t_config)) {
    auto event = sdsl::memory_monitor::event("Successor");
    const auto n = sdsl::int_vector_buffer<>(cache_file_name(keys[kBWT][kBase], t_config)).size();
    constructBitVectorFromIntVector<typename Index::BvMarks>(keys[kPsi][kTail][kTextPos], t_config, n, false, true);
  }

  // Construct Links from Mark to Sample
  if (!sdsl::cache_file_exists<typename Index::MarksToSamples>(keys[kPsi][kTail][kTextPosAsc][kLink], t_config)) {
    auto event = sdsl::memory_monitor::event("Mark2Sample Links");

    sdsl::int_vector<> mark_to_sample_links;
    if (!sdsl::cache_file_exists<sdsl::int_vector<>>(keys[kPsi][kTail][kTextPosAsc][kLink], t_config)) {
      mark_to_sample_links = constructMarkToSampleLinksForPhiForwardWithPsiRuns(t_config);
    } else {
      sdsl::load_from_cache(mark_to_sample_links, keys[kPsi][kTail][kTextPosAsc][kLink], t_config, true);
    }

    if (!std::is_same_v<typename Index::MarksToSamples, sdsl::int_vector<>>) {
      auto values = sri::construct<typename Index::MarksToSamples>(mark_to_sample_links);
      sri::store_to_cache(values, keys[kPsi][kTail][kTextPosAsc][kLink], t_config, true);
    }
  }
}

}  // namespace sri

#endif  // SRI_R_CSA_H_
