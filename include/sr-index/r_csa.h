//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/15/21.
//

#ifndef SRI_CSA_H_
#define SRI_CSA_H_

#include <map>
#include <string>
#include <any>
#include <functional>
#include <memory>

#include <sdsl/csa_alphabet_strategy.hpp>

#include "alphabet.h"
#include "psi.h"
#include "lf.h"
#include "construct.h"
#include "sequence_ops.h"
#include "index_base.h"

namespace sri {

template<typename TStorage = GenericStorage,
    typename TAlphabet = Alphabet<>,
    typename TPsiRLE = PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>>
class RCSA : public IndexBaseWithExternalStorage<TStorage> {
 public:
  using Base = IndexBaseWithExternalStorage<TStorage>;

  explicit RCSA(const TStorage &t_storage) : Base(t_storage) {}

  RCSA() = default;

  void load(sdsl::cache_config t_config) override {
    TSource source(std::ref(t_config));
    loadInner(source);
  }

  void load(std::istream &in) override {
    TSource source(std::ref(in));
    loadInner(source);
  }

  using typename Base::size_type;
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(key(SrIndexKey::ALPHABET), out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), out, child, "psi");

    written_bytes += this->template serializeItem<TSample>(key(SrIndexKey::SAMPLES), out, child, "samples");

    written_bytes += this->template serializeItem<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_select");

    written_bytes += this->template serializeItem<TMarkToSampleIdx>(
        key(SrIndexKey::MARK_TO_SAMPLE), out, child, "mark_to_sample");

    return written_bytes;
  }

 protected:

  using typename Base::TSource;

  virtual void loadInner(TSource &t_source) {
    setupKeyNames();
    loadAllItems(t_source);
    constructIndex(t_source);
  }

  using Base::key;
  virtual void setupKeyNames() {
    if (!this->keys_.empty()) return;

    this->keys_.resize(5);
    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
    key(SrIndexKey::SAMPLES) = conf::KEY_BWT_RUN_FIRST_TEXT_POS;
    key(SrIndexKey::MARKS) = conf::KEY_BWT_RUN_LAST_TEXT_POS;
    key(SrIndexKey::MARK_TO_SAMPLE) = conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
  }

  virtual void loadAllItems(TSource &t_source) {
    this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);

    this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);

    this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);

    this->template loadItem<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    this->template loadBVRank<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    this->template loadBVSelect<TBvMark>(key(SrIndexKey::MARKS), t_source, true);

    this->template loadItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), t_source);
  }

  virtual void constructIndex(TSource &t_source) {
    this->index_.reset(new RIndexBase{
        constructLF(t_source),
        constructComputeDataBackwardSearchStep(
            [](const Range &tt_range, auto tt_c, const RangeLF &tt_next_range, std::size_t tt_step) {
              const auto &[start, end] = tt_next_range;
              return DataBackwardSearchStep{tt_step, RunData(start.run.start)};
            }),
        constructComputeSAValues(constructPhiForRange(t_source), constructComputeToehold(t_source)),
        this->n_,
        [](const auto &tt_step) { return DataBackwardSearchStep{0, RunData(0)}; },
        constructGetSymbol(t_source),
        [](auto tt_seq_size) { return Range{0, tt_seq_size}; },
        constructIsRangeEmpty()
    });
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

  using Char = typename TAlphabet::comp_char_type;
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

  auto constructIsRangeEmpty() {
    return [](const Range &tt_range) { return !(tt_range.start < tt_range.end); };
  }

  auto constructLF(TSource &t_source) {
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);
    auto cumulative = RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));
    this->n_ = cumulative[cref_alphabet.get().sigma];

    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
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

    return LF(psi_rank, cumulative, create_range, empty_range);
  }

  struct RunData {
    std::size_t pos;

    explicit RunData(std::size_t t_pos) : pos{t_pos} {}
    virtual ~RunData() = default;
  };

  using DataBackwardSearchStep = sri::DataBackwardSearchStep<RunData>;

  template<typename TCreateData>
  auto constructComputeDataBackwardSearchStep(const TCreateData &t_create_data) {
    return buildComputeDataBackwardSearchStepForPhiForward(t_create_data);
  }

  using Value = std::size_t;
  auto constructPhiForRange(TSource &t_source) {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto successor = CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), t_source);
    auto get_mark_to_sample_idx = RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, true);

    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto get_sample = RandomAccessForCRefContainer(cref_samples);
    SampleValidatorDefault sample_validator_default;

    auto phi = buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);
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

  auto constructComputeToehold(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);

    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto get_sa_value_for_bwt_run_start = [cref_psi_core, cref_samples](const RunData &tt_run_data) {
      auto run = cref_psi_core.get().rankRun(tt_run_data.pos);
      return cref_samples.get()[run] + 1;
    };

    return buildComputeToeholdForPhiForward(get_sa_value_for_bwt_run_start, cref_psi_core.get().size());
  }

  template<typename TPhiRange, typename TComputeToehold>
  auto constructComputeSAValues(const TPhiRange &t_phi_range, const TComputeToehold &t_compute_toehold) {
    auto update_range = [](Range tt_range) {
      auto &[start, end] = tt_range;
      ++start;
      return tt_range;
    };

    return ComputeAllValuesWithPhiForRange(t_phi_range, t_compute_toehold, update_range);
  }

  auto constructGetSymbol(TSource &t_source) {
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);

    auto get_symbol = [cref_alphabet](char tt_c) { return cref_alphabet.get().char2comp[(uint8_t) tt_c]; };
    return get_symbol;
  }

};

template<typename TStorage = GenericStorage,
    typename TAlphabet = Alphabet<>,
    typename TPsiRLE = PsiCoreRLE<>,
    typename TSA = sdsl::int_vector<>>
class CSARaw : public RCSA<TStorage, TAlphabet, TPsiRLE> {
 public:
  using Base = RCSA<TStorage, TAlphabet, TPsiRLE>;

  explicit CSARaw(const TStorage &t_storage) : Base(t_storage) {}

  CSARaw() = default;

  using typename Base::size_type;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(conf::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += this->template serializeItem<TSA>(sdsl::conf::KEY_SA, out, child, "sa");

    return written_bytes;
  }

 protected:

  using typename Base::TSource;
  using Base::key;

  virtual void loadAllItems(TSource &t_source) {
    this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);
    this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    this->template loadItem<TSA>(sdsl::conf::KEY_SA, t_source);
  }

  using typename Base::Range;
  using typename Base::RangeLF;
  using typename Base::DataBackwardSearchStep;
  using typename Base::RunData;
  virtual void constructIndex(TSource &t_source) {
    this->index_.reset(new RIndexBase{
        this->constructLF(t_source),
        this->constructComputeDataBackwardSearchStep(
            [](const Range &tt_range, auto tt_c, const RangeLF &tt_next_range, std::size_t tt_step) {
              const auto &[start, end] = tt_next_range;
              return DataBackwardSearchStep{tt_step, RunData(start.run.start)};
            }),
        constructComputeSAValues(t_source),
        this->n_,
        [](const auto &tt_step) { return DataBackwardSearchStep{0, RunData(0)}; },
        this->constructGetSymbol(t_source),
        [](auto tt_seq_size) { return Range{0, tt_seq_size}; },
        this->constructIsRangeEmpty()
    });
  }

  auto constructComputeSAValues(TSource &t_source) {
    auto cref_sa = this->template loadItem<TSA>(sdsl::conf::KEY_SA, t_source);

    auto compute_sa_values = [cref_sa](const auto &tt_range, const auto &, auto tt_report) {
      auto [first, last] = tt_range;

      while (first < last) {
        tt_report(cref_sa.get()[first++]);
      }
    };

    return compute_sa_values;
  }

};

template<uint8_t t_width, typename TBvMark>
void constructCSA(const std::string &t_data_path, sdsl::cache_config &t_config);

template<typename TStorage, template<uint8_t> typename TAlphabet, uint8_t t_width, typename TPsiCore, typename TBvMark, typename TMarkToSampleIdx, typename TSample>
void construct(RCSA<TStorage, TAlphabet<t_width>, TPsiCore, TBvMark, TMarkToSampleIdx, TSample> &t_index,
               const std::string &t_data_path,
               sdsl::cache_config &t_config) {
  constructCSA<t_width, TBvMark>(t_data_path, t_config);

  t_index.load(t_config);
}

template<uint8_t t_width, typename TBvMark>
void constructCSA(const std::string &t_data_path, sdsl::cache_config &t_config) {
  constructIndexBaseItems<t_width>(t_data_path, t_config);

  {
    // Construct Psi
    auto event = sdsl::memory_monitor::event("Psi");
    if (!cache_file_exists(sdsl::conf::KEY_PSI, t_config)) {
      constructPsi<t_width>(t_config);
    }
  }

  {
    // Construct Links from Mark to Sample
    auto event = sdsl::memory_monitor::event("Mark2Sample Links");
    if (!cache_file_exists(conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config)) {
      constructMarkToSampleLinksForPhiForward<t_width>(t_config);
    }
  }

  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = bwt_buf.size();
  }

  {
    // Construct Successor on the text positions of BWT run last letter
    auto event = sdsl::memory_monitor::event("Successor");
    if (!sdsl::cache_file_exists<TBvMark>(conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config)) {
      constructBitVectorFromIntVector<TBvMark>(conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config, n, false);
    }
  }
}

}

#endif //SRI_CSA_H_