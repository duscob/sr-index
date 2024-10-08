//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_SR_CSA_H_
#define SRI_SR_CSA_H_

#include <optional>

#include "sampling.h"
#include "r_csa.h"

namespace sri {

template<typename TStorage = GenericStorage,
    typename TAlphabet = Alphabet<>,
    typename TPsiRLE = PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>>
class SrCSABase : public RCSAWithBWTRun<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample> {
 public:
  using Base = RCSAWithBWTRun<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample>;

  SrCSABase(const TStorage &t_storage, std::size_t t_sr)
      : Base(t_storage), subsample_rate_{t_sr}, key_prefix_{std::to_string(subsample_rate_) + "_"} {
  }

  explicit SrCSABase(std::size_t t_sr)
      : Base(), subsample_rate_{t_sr}, key_prefix_{std::to_string(subsample_rate_) + "_"} {}

  SrCSABase() = default;

  virtual ~SrCSABase() = default;

  [[nodiscard]] std::size_t SubsampleRate() const { return subsample_rate_; }

  void load(Config t_config) override {
    TSource source(std::ref(t_config));
    this->loadInner(source);
  }

  void load(std::istream &in) override {
    sdsl::read_member(this->subsample_rate_, in);
    TSource source(std::ref(in));
    this->loadInner(source);
  }

  using typename Base::size_type;
  using typename Base::ItemKey;
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += sdsl::write_member(this->subsample_rate_, out, child, "subsample_rate");

    written_bytes += Base::serialize(out, v, name);

    return written_bytes;
  }

 protected:

  using Base::key;
  void setupKeyNames() override {
    Base::setupKeyNames();
//    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
    key(ItemKey::SAMPLES) = key_prefix_ + conf::KEY_BWT_RUN_FIRST_TEXT_POS;
    key(ItemKey::MARKS) = key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
    key(ItemKey::MARK_TO_SAMPLE) = key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
    key(ItemKey::SAMPLES_IDX) = key_prefix_ + conf::KEY_BWT_RUN_FIRST_IDX;
  }

  using typename Base::TSource;

  using typename Base::Char;

  template<typename TCreateRun>
  auto constructSplitRunInBWTRuns(TSource &t_source, const TCreateRun &t_create_run) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);

    return [cref_psi_core, cref_alphabet, t_create_run](const auto &tt_run) -> auto {
      const auto &first = tt_run.start;
      const auto &last = tt_run.end;
      const auto &cumulative = cref_alphabet.get().C;
      auto c = computeCForSAIndex(cumulative, first);

      auto cum_c = cumulative[c];
      auto cum_next_c = cumulative[c + 1];
      auto first_rank = first - cum_c + 1;
      auto last_rank = std::min(last, cum_next_c) - cum_c + 1;
      std::vector<decltype(t_create_run(0u, 0u, (Char) 0u, 0u, false))> runs;
      auto report = [&runs, &t_create_run](auto tt_first, auto tt_last, auto tt_c, auto tt_n_run, auto tt_is_first) {
        runs.emplace_back(t_create_run(tt_first, tt_last, tt_c, tt_n_run, tt_is_first));
      };
      cref_psi_core.get().computeForwardRuns(c, first_rank, last_rank, report);

      while (cum_next_c < last) {
        ++c;
        cum_c = cum_next_c;
        cum_next_c = cumulative[c + 1];
        last_rank = std::min(last, cum_next_c) - cum_c + 1;
        cref_psi_core.get().computeForwardRuns(c, 1, last_rank, report);
      }

      return runs;
    };
  }

  using BvMark = TBvMark;
  using MarkToSampleIdx = TMarkToSampleIdx;
  using Sample = TSample;

  template<typename TGetSampleRun, typename TSplitRangeInBWTRuns, typename TSplitRunInBWTRuns, typename TUpdateRun, typename TIsRunEmpty>
  auto constructPhiForRange(TSource &t_source,
                            const TGetSampleRun &t_get_sample,
                            const TSplitRangeInBWTRuns &t_split_range,
                            const TSplitRunInBWTRuns &t_split_run,
                            const TUpdateRun &t_update_run,
                            const TIsRunEmpty &t_is_run_empty) {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(key(ItemKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(key(ItemKey::MARKS), t_source, true);
    auto successor = CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(key(ItemKey::MARK_TO_SAMPLE), t_source);
    auto get_mark_to_sample_idx = RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, false);

    auto cref_samples = this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source);
    auto get_sample = RandomAccessForCRefContainer(cref_samples);
    SampleValidatorDefault sample_validator_default;

    auto phi = buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);
    auto phi_simple = [phi](const auto &tt_prev_value) { return phi(tt_prev_value).first; };

    return PhiForwardForRange(phi_simple,
                              t_get_sample,
                              t_split_range,
                              t_split_run,
                              subsample_rate_,
                              this->n_,
                              this->constructIsRangeEmpty(),
                              t_update_run,
                              t_is_run_empty);
  }

  using typename Base::RunData;

  template<typename TGetSampleRunData, typename TPsiRunData>
  auto constructComputeToehold(TSource &t_source, const TGetSampleRunData &t_get_sample, const TPsiRunData &t_psi) {
    auto compute_sa_value_run_data = buildComputeSAValueForward(t_get_sample, t_psi, this->n_);
    return buildComputeToeholdForPhiForward(compute_sa_value_run_data, this->n_);
  }

  std::size_t subsample_rate_ = 1;
  std::string key_prefix_;
};

template<typename TStorage = GenericStorage,
    typename TAlphabet = Alphabet<>,
    typename TPsiRLE = PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSampleIdx = sdsl::sd_vector<>,
    typename TRunCumulativeCount = sdsl::int_vector<>>
class SrCSASlim : public SrCSABase<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample> {
 public:
  using Base = SrCSABase<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample>;

  SrCSASlim(const TStorage &t_storage, std::size_t t_sr) : Base(t_storage, t_sr) {}

  explicit SrCSASlim(std::size_t t_sr) : Base(t_sr) {}

  SrCSASlim() = default;

  using typename Base::size_type;
  using typename Base::ItemKey;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += Base::serialize(out, v, name);

    written_bytes += this->template serializeItem<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), out, child, "samples_idx");
    written_bytes +=
        this->template serializeRank<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), out, child, "samples_idx_rank");

    written_bytes += this->template serializeItem<TRunCumulativeCount>(key(ItemKey::RUN_CUMULATIVE_COUNT),
                                                                       out,
                                                                       child,
                                                                       "run_cumulative_count");

    return written_bytes;
  }

 protected:

  using Base::key;

  void setupKeyNames() override {
    Base::setupKeyNames();
//    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
    key(ItemKey::SAMPLES) = KeySortedByAlphabet(this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_TEXT_POS);
    key(ItemKey::MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
    key(ItemKey::MARK_TO_SAMPLE) = KeySortedByAlphabet(
        this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX);
    key(ItemKey::SAMPLES_IDX) = KeySortedByAlphabet(this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_IDX);
    key(ItemKey::RUN_CUMULATIVE_COUNT) = conf::KEY_BWT_RUN_CUMULATIVE_COUNT;
  }

  using typename Base::TSource;

  void loadAllItems(TSource &t_source) override {
    Base::loadAllItems(t_source);

    this->template loadItem<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), t_source, true);
    this->template loadBVRank<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), t_source, true);

    this->template loadItem<TRunCumulativeCount>(key(ItemKey::RUN_CUMULATIVE_COUNT), t_source);
  }

  void constructIndex(TSource &t_source) override {
    constructIndex(t_source,
                   [this](auto &t_source) {
                     return this->constructPhiForRange(t_source,
                                                       constructGetSampleForRun(t_source),
                                                       constructSplitRangeInBWTRuns(t_source),
                                                       constructSplitRunInBWTRuns(t_source),
                                                       constructUpdateRun(),
                                                       constructIsRunEmpty());
                   });
  }

  using typename Base::Range;
  using typename Base::RangeLF;
  template<typename TPhiRange>
  void constructIndex(TSource &t_source, const TPhiRange &t_phi_range) {
    this->index_.reset(new RIndexBase{
        this->constructLF(t_source),
        this->constructComputeDataBackwardSearchStep(
            [](const Range &tt_range, Char tt_c, const RangeLF &tt_next_range, std::size_t tt_step) {
              const auto &run = tt_next_range.start.run;
              return DataBackwardSearchStep{tt_step,
                                            RunDataExt(run.start, tt_c, run.rank, tt_range.start != run.start)};
            }),
        this->constructComputeSAValues(
            t_phi_range(t_source),
            this->constructComputeToehold(t_source,
                                          constructGetSampleForRunData(t_source),
                                          constructPsiForRunData(t_source))),
        this->n_,
        [](const auto &tt_step) { return DataBackwardSearchStep{0, RunDataExt()}; },
        this->constructGetSymbol(t_source),
        [](auto tt_seq_size) { return Range{0, tt_seq_size}; },
        this->constructIsRangeEmpty()
    });
  }

  using typename Base::RunData;
  using typename Base::Char;
  struct RunDataExt : RunData {
    Char c;
    std::size_t partial_rank;
    bool is_run_start;

    explicit RunDataExt(std::size_t t_pos = 0, Char t_c = 0, std::size_t t_partial_rank = 0, bool t_is_run_start = false)
        : RunData(t_pos), c{t_c}, partial_rank{t_partial_rank}, is_run_start{t_is_run_start} {}
  };

  using DataBackwardSearchStep = sri::DataBackwardSearchStep<RunDataExt>;

  using typename Base::Value;
  auto constructGetSample(TSource &t_source) {
    auto cref_run_cum_c = this->template loadItem<TRunCumulativeCount>(key(ItemKey::RUN_CUMULATIVE_COUNT), t_source);
    auto cref_bv_sample_idx = this->template loadItem<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_idx_rank = this->template loadBVRank<TBvSampleIdx>(key(ItemKey::SAMPLES_IDX), t_source, true);
    auto cref_samples = this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source);

    auto get_sample = [cref_run_cum_c, cref_bv_sample_idx, bv_sample_idx_rank, cref_samples](
        auto tt_char, auto tt_partial_rank, bool tt_is_run_start) -> std::optional<Value> {
      if (!tt_is_run_start) return std::nullopt;

      std::size_t cc = 0 < tt_char ? cref_run_cum_c.get()[tt_char - 1] : 0;
      auto idx = cc + tt_partial_rank;
      if (!cref_bv_sample_idx.get()[idx]) return std::nullopt;

      auto rank = bv_sample_idx_rank(idx);
      return cref_samples.get()[rank];
    };

    return get_sample;
  }

  auto constructGetSampleForRunData(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);
    return [get_sample](const RunDataExt &tt_run_data) {
      return get_sample(tt_run_data.c, tt_run_data.partial_rank, tt_run_data.is_run_start);
    };
  }

  struct Run {
    std::size_t start = 0;
    std::size_t end = 0;
    Char c = 0;
    std::size_t partial_rank = 0;
    bool is_run_start = false;

    bool operator<(const Run &rhs) const {
      return start < rhs.start;
    }
  };

  auto constructUpdateRun() {
    return [](Run &tt_run) {
      ++tt_run.start;
      tt_run.is_run_start = false;
    };
  }

  auto constructIsRunEmpty() {
    return [](const Run &tt_run) { return !(tt_run.start < tt_run.end); };
  }

  auto constructGetSampleForRun(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);

    return [get_sample](const Run &tt_run) {
      return get_sample(tt_run.c, tt_run.partial_rank, tt_run.is_run_start);
    };
  }

  auto constructCreateRun() {
    auto create_run = [](auto tt_first, auto tt_last, auto tt_c, auto tt_partial_rank, auto tt_is_first) {
      return Run{tt_first, tt_last, tt_c, tt_partial_rank, tt_is_first};
    };

    return create_run;
  }

  auto constructSplitRangeInBWTRuns(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);
    auto create_run = constructCreateRun();

    return [cref_psi_core, create_run](const auto &tt_range) -> auto {
      const auto &[first, last] = tt_range;
      return cref_psi_core.get().splitInSortedRuns(first, last, create_run);
    };
  }

  auto constructSplitRunInBWTRuns(TSource &t_source) {
    return Base::constructSplitRunInBWTRuns(t_source, constructCreateRun());
  }

  auto constructPsiForRunData(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);
    auto psi_select = [cref_psi_core](Char tt_c, auto tt_rnk) {
      RunDataExt run_data;
      auto report = [&run_data, &tt_c](auto tt_pos, auto tt_run_start, auto tt_run_end, auto tt_run_rank) {
        run_data = RunDataExt{tt_pos, tt_c, tt_run_rank, tt_pos == tt_run_start};
      };
      cref_psi_core.get().select(tt_c, tt_rnk, report);
      return run_data;
    };

    auto cref_alphabet = this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);
    auto get_c = [cref_alphabet](auto tt_index) { return computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));

    auto psi = Psi(psi_select, get_c, cumulative);

    return [psi](const RunDataExt &tt_run_data) { return psi(tt_run_data.pos); };
  }
};

template<typename TStorage = GenericStorage,
    typename TAlphabet = Alphabet<>,
    typename TPsiRLE = PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSamplePos = sdsl::sd_vector<>>
class SrCSA : public SrCSABase<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample> {
 public:
  using Base = SrCSABase<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample>;

  SrCSA(const TStorage &t_storage, std::size_t t_sr) : Base(t_storage, t_sr) {}

  explicit SrCSA(std::size_t t_sr) : Base(t_sr) {}

  SrCSA() = default;

  using typename Base::size_type;
  using typename Base::ItemKey;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += Base::serialize(out, v, name);

    written_bytes += this->template serializeItem<TBvSamplePos>(
        key(ItemKey::SAMPLES_IDX), out, child, "samples_idx");
    written_bytes += this->template serializeRank<TBvSamplePos>(
        key(ItemKey::SAMPLES_IDX), out, child, "samples_idx_rank");
    written_bytes += this->template serializeSelect<TBvSamplePos>(
        key(ItemKey::SAMPLES_IDX), out, child, "samples_idx_select");

    return written_bytes;
  }

 protected:

  using Base::key;
  void setupKeyNames() override {
    Base::setupKeyNames();
//    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
//    key(SrIndexKey::SAMPLES) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_TEXT_POS;
//    key(SrIndexKey::MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
//    key(SrIndexKey::MARK_TO_SAMPLE) =
//        this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
    key(ItemKey::SAMPLES_IDX) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST;
  }

  using typename Base::TSource;

  void loadAllItems(TSource &t_source) override {
    Base::loadAllItems(t_source);

    this->template loadItem<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);
    this->template loadBVRank<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);
  }

  void constructIndex(TSource &t_source) override {
    constructIndex(t_source,
                   [this](auto &t_source) {
                     return this->constructPhiForRange(t_source,
                                                       constructGetSampleForRun(t_source),
                                                       constructSplitRangeInBWTRuns(t_source),
                                                       constructSplitRunInBWTRuns(t_source),
                                                       constructUpdateRun(),
                                                       constructIsRunEmpty());
                   });
  }

  using typename Base::Range;
  using typename Base::RangeLF;
  using typename Base::DataBackwardSearchStep;
  using typename Base::RunData;
  template<typename TPhiRange>
  void constructIndex(TSource &t_source, const TPhiRange &t_phi_range) {
    this->index_.reset(new RIndexBase{
        this->constructLF(t_source),
        this->constructComputeDataBackwardSearchStep(
            [](const auto &tt_range, auto tt_c, const RangeLF &tt_next_range, std::size_t tt_step) {
              const auto &[start, end] = tt_next_range;
              return DataBackwardSearchStep{tt_step, RunData(start.run.start)};
            }),
        this->constructComputeSAValues(
            t_phi_range(t_source),
            this->constructComputeToehold(t_source,
                                          constructGetSampleForRunData(t_source),
                                          constructPsiForRunData(t_source))),
        this->n_,
        [](const auto &tt_step) { return DataBackwardSearchStep{0, RunData(0)}; },
        this->constructGetSymbol(t_source),
        [](auto tt_seq_size) { return Range{0, tt_seq_size}; },
        this->constructIsRangeEmpty()
    });
  }

  auto constructGetSample(TSource &t_source) {
    auto cref_samples = this->template loadItem<TSample>(key(ItemKey::SAMPLES), t_source);
    auto cref_bv_sample_pos = this->template loadItem<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);

    auto get_sample = [cref_samples, cref_bv_sample_pos, bv_sample_pos_rank](auto tt_pos) {
      return cref_bv_sample_pos.get()[tt_pos]
             ? std::optional<std::size_t>(cref_samples.get()[bv_sample_pos_rank(tt_pos)])
             : std::nullopt;
    };

    return get_sample;
  }

  auto constructGetSampleForRunData(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);
    return [get_sample](const RunData &tt_run_data) { return get_sample(tt_run_data.pos); };
  }

  auto constructGetSampleForRun(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);

    return [get_sample](const auto &tt_run) {
      return get_sample(tt_run.start);
    };
  }

  struct Run {
    std::size_t start = 0;
    std::size_t end = 0;
  };

  auto constructUpdateRun() {
    return [](Run &tt_run) { ++tt_run.start; };
  }

  auto constructIsRunEmpty() {
    return [](const Run &tt_run) { return !(tt_run.start < tt_run.end); };
  }

  auto constructCreateRun() {
    return [](auto tt_first, auto tt_last, auto, auto, auto) { return Run{tt_first, tt_last}; };
  }

  auto constructSplitRangeInBWTRuns(TSource &t_source) {
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_pos_select = this->template loadBVSelect<TBvSamplePos>(key(ItemKey::SAMPLES_IDX), t_source, true);

    return [bv_sample_pos_rank, bv_sample_pos_select](const auto &tt_range) -> auto {
      const auto &[first, last] = tt_range;
      auto first_rank = bv_sample_pos_rank(first + 1) + 1;
      auto last_rank = bv_sample_pos_rank(last);

      std::vector<Run> runs;
      auto run_start = first;
      for (int i = first_rank; i <= last_rank; ++i) {
        auto run_end = bv_sample_pos_select(i);
        runs.emplace_back(Run{run_start, run_end});

        run_start = run_end;
      }

      runs.emplace_back(Run{run_start, last});

      return runs;
    };
  }

  auto constructSplitRunInBWTRuns(TSource &t_source) {
    return Base::constructSplitRunInBWTRuns(t_source, constructCreateRun());
  }

  auto constructPsiForRunData(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(ItemKey::NAVIGATE), t_source, true);
    auto psi_select = [cref_psi_core](auto tt_c, auto tt_rnk) { return cref_psi_core.get().select(tt_c, tt_rnk); };

    auto cref_alphabet = this->template loadItem<TAlphabet>(key(ItemKey::ALPHABET), t_source);
    auto get_c = [cref_alphabet](auto tt_index) { return computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));

    auto psi = Psi(psi_select, get_c, cumulative);

    return [psi](RunData tt_run_data) {
      tt_run_data.pos = psi(tt_run_data.pos);
      return tt_run_data;
    };
  }
};

template<typename TSrCSA, typename TBvValidMark = sdsl::bit_vector>
class SrCSAValidMark : public TSrCSA {
 public:
  using Base = TSrCSA;

  template<typename TStorage>
  SrCSAValidMark(const TStorage &t_storage, std::size_t t_sr) : Base(t_storage, t_sr) {}

  explicit SrCSAValidMark(std::size_t t_sr) : Base(t_sr) {}

  SrCSAValidMark() = default;

  using typename Base::size_type;
  using typename Base::ItemKey;
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = Base::serialize(out, v, name);

    written_bytes += this->template serializeItem<TBvValidMark>(
        key(ItemKey::VALID_MARKS), out, child, "valid_marks");

    return written_bytes;
  }

 protected:

  using Base::key;

  void setupKeyNames() override {
    Base::setupKeyNames();
//    constexpr uint8_t t_width = SrCSAValidMark<TSrCSA>::AlphabetWidth;
//    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
//    key(SrIndexKey::SAMPLES) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_TEXT_POS;
//    key(SrIndexKey::MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
//    key(SrIndexKey::MARK_TO_SAMPLE) =
//        this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
//    key(SrIndexKey::SAMPLES_IDX) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_IDX;
    key(ItemKey::VALID_MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK;
  }

  using typename Base::TSource;

  void loadAllItems(TSource &t_source) override {
    Base::loadAllItems(t_source);

    this->template loadItem<TBvValidMark>(key(ItemKey::VALID_MARKS), t_source, true);
  }

  void constructIndex(TSource &t_source) override {
    Base::constructIndex(t_source,
                         [this](auto &t_source) {
                           return constructPhiForRange(t_source,
                                                       Base::constructGetSampleForRun(t_source),
                                                       Base::constructSplitRangeInBWTRuns(t_source),
                                                       Base::constructSplitRunInBWTRuns(t_source),
                                                       Base::constructUpdateRun(),
                                                       Base::constructIsRunEmpty(),
                                                       []() { return SampleValidatorDefault(); });
                         });
  }

  using Base::constructIndex;

  using typename Base::BvMark;
  using typename Base::MarkToSampleIdx;
  using typename Base::Sample;
  template<typename TGetSampleRun, typename TSplitRangeInBWTRuns, typename TSplitRunInBWTRuns, typename TUpdateRun, typename TIsRunEmpty, typename TConstructValidateSample>
  auto constructPhiForRange(TSource &t_source,
                            const TGetSampleRun &t_get_sample,
                            const TSplitRangeInBWTRuns &t_split_range,
                            const TSplitRunInBWTRuns &t_split_run,
                            const TUpdateRun &t_update_run,
                            const TIsRunEmpty &t_is_run_empty,
                            const TConstructValidateSample &t_construct_validate_sample) {
    auto bv_mark_rank = this->template loadBVRank<BvMark>(key(ItemKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<BvMark>(key(ItemKey::MARKS), t_source, true);
    auto successor = CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<MarkToSampleIdx>(key(ItemKey::MARK_TO_SAMPLE), t_source);
    auto cref_bv_valid_mark = this->template loadItem<TBvValidMark>(key(ItemKey::VALID_MARKS), t_source, true);
    auto get_mark_to_sample_idx = RandomAccessForTwoContainers(cref_mark_to_sample_idx, cref_bv_valid_mark);

    auto cref_samples = this->template loadItem<Sample>(key(ItemKey::SAMPLES), t_source);
    auto get_sample = RandomAccessForCRefContainer(cref_samples);

    auto phi = buildPhiForward(successor, get_mark_to_sample_idx, get_sample, t_construct_validate_sample(), this->n_);

    return PhiForwardForRangeWithValidity(phi,
                                          t_get_sample,
                                          t_split_range,
                                          t_split_run,
                                          this->subsample_rate_,
                                          this->n_,
                                          this->constructIsRangeEmpty(),
                                          t_update_run,
                                          t_is_run_empty);
  }
};

template<typename TSrCSA, typename TBvValidMark = sdsl::bit_vector, typename TValidArea = sdsl::int_vector<>>
class SrCSAValidArea : public SrCSAValidMark<TSrCSA, TBvValidMark> {
 public:
  using Base = SrCSAValidMark<TSrCSA, TBvValidMark>;

  template<typename TStorage>
  SrCSAValidArea(const TStorage &t_storage, std::size_t t_sr) : Base(t_storage, t_sr) {}

  explicit SrCSAValidArea(std::size_t t_sr) : Base(t_sr) {}

  SrCSAValidArea() = default;

  using typename Base::size_type;
  using typename Base::ItemKey;
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = Base::serialize(out, v, name);

    written_bytes += this->template serializeRank<TBvValidMark, typename TBvValidMark::rank_0_type>(
        key(ItemKey::VALID_MARKS), out, child, "valid_marks_rank");

    written_bytes += this->template serializeItem<TValidArea>(key(ItemKey::VALID_AREAS), out, child, "valid_areas");

    return written_bytes;
  }

 protected:

  using Base::key;
  void setupKeyNames() override {
    Base::setupKeyNames();
//    constexpr uint8_t t_width = SrCSAValidMark<TSrCSA>::AlphabetWidth;
//    key(SrIndexKey::ALPHABET) = conf::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
//    key(SrIndexKey::SAMPLES) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_TEXT_POS;
//    key(SrIndexKey::MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
//    key(SrIndexKey::MARK_TO_SAMPLE) =
//        this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
//    key(SrIndexKey::SAMPLES_IDX) = this->key_prefix_ + conf::KEY_BWT_RUN_FIRST_IDX;
//    key(SrIndexKey::VALID_MARKS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK;
    key(ItemKey::VALID_AREAS) = this->key_prefix_ + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA;
  }

  using typename Base::TSource;

  void loadAllItems(TSource &t_source) override {
    Base::loadAllItems(t_source);

    this->template loadBVRank<TBvValidMark, typename TBvValidMark::rank_0_type>(
        key(ItemKey::VALID_MARKS), t_source, true);

    this->template loadItem<TValidArea>(key(ItemKey::VALID_AREAS), t_source);
  }

  void constructIndex(TSource &t_source) override {
    Base::constructIndex(t_source,
                         [this](auto &t_source) {
                           return Base::constructPhiForRange(t_source,
                                                             Base::constructGetSampleForRun(t_source),
                                                             Base::constructSplitRangeInBWTRuns(t_source),
                                                             Base::constructSplitRunInBWTRuns(t_source),
                                                             Base::constructUpdateRun(),
                                                             Base::constructIsRunEmpty(),
                                                             constructValidateSample(t_source));
                         });
  }

  auto constructValidateSample(TSource &t_source) {
    return [&t_source, this]() {
      auto bv_valid_mark_rank = this->template loadBVRank<TBvValidMark, typename TBvValidMark::rank_0_type>(
          key(ItemKey::VALID_MARKS), t_source, true);
      auto cref_valid_area = this->template loadItem<TValidArea>(key(ItemKey::VALID_AREAS), t_source);
      auto get_valid_area = RandomAccessForCRefContainer(cref_valid_area);

      return SampleValidator(bv_valid_mark_rank, get_valid_area);
    };
  }
};

template<uint8_t t_width, typename TBVMark>
void constructSrCSACommonsWithBWTRuns(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

void constructSamplesSortedByAlphabet(sdsl::cache_config &t_config);

void constructSubsamplingBackwardSamplesSortedByAlphabet(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<typename TStorage, template<uint8_t> typename TAlphabet, uint8_t t_width, typename TPsiCore, typename TBvMark,
  typename TMarkToSample, typename TSample, typename TBVSampleIdx, typename TRunCumCnt>
void construct(
  SrCSASlim<TStorage, TAlphabet<t_width>, TPsiCore, TBvMark, TMarkToSample, TSample, TBVSampleIdx, TRunCumCnt>& t_index,
  const std::string& t_data_path,
  Config& t_config
) {
  constructRCSAWithBWTRuns<t_width, TBvMark>(t_data_path, t_config);

  auto subsample_rate = t_index.SubsampleRate();

  constructSrCSACommonsWithBWTRuns<t_width, TBvMark>(subsample_rate, t_config);

  // Construct samples' indices sorted by alphabet
  if (!sdsl::cache_file_exists(KeySortedByAlphabet(conf::KEY_BWT_RUN_FIRST_IDX), t_config)) {
    auto event = sdsl::memory_monitor::event("Samples");
    constructSamplesSortedByAlphabet(t_config);
  }

  auto prefix_key = std::to_string(subsample_rate) + "_";

  // Construct subsampling backward of samples sorted by alphabet
  if (!sdsl::cache_file_exists(KeySortedByAlphabet(prefix_key + conf::KEY_BWT_RUN_FIRST_TEXT_POS), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardSamplesSortedByAlphabet(subsample_rate, t_config);
  }

  // Construct subsampling indices backward of samples sorted by alphabet
  if (!sdsl::cache_file_exists<TBVSampleIdx>(KeySortedByAlphabet(prefix_key + conf::KEY_BWT_RUN_FIRST_IDX), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    const auto r = sdsl::int_vector_buffer<>(sdsl::cache_file_name(conf::KEY_BWT_RUN_FIRST, t_config)).size();

    constructBitVectorFromIntVector<TBVSampleIdx>(KeySortedByAlphabet(prefix_key + conf::KEY_BWT_RUN_FIRST_IDX),
                                                  t_config,
                                                  r,
                                                  false);
  }

  t_index.load(t_config);
}

void constructSubsamplingBackwardSamplesPosition(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<typename TStorage, template<uint8_t> typename TAlphabet, uint8_t t_width, typename TPsiCore, typename TBvMark, typename TMarkToSampleIdx, typename TSample, typename TBvSamplePos>
void construct(SrCSA<TStorage, TAlphabet<t_width>, TPsiCore, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos> &t_index,
               const std::string &t_data_path,
               sri::Config &t_config) {
  constructRCSAWithBWTRuns<t_width, TBvMark>(t_data_path, t_config);

  auto subsample_rate = t_index.SubsampleRate();

  constructSrCSACommonsWithBWTRuns<t_width, TBvMark>(subsample_rate, t_config);

  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = buf.size();
  }

  auto prefix = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling backward of samples (text positions of BWT-run last letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = prefix + conf::KEY_BWT_RUN_FIRST;
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSubsamplingBackwardSamplesPosition(subsample_rate, t_config);
    }

    if (!sdsl::cache_file_exists<TBvSamplePos>(key, t_config)) {
      constructBitVectorFromIntVector<TBvSamplePos>(key, t_config, n, false);
    }
  }

  t_index.load(t_config);
}

void constructSubsamplingBackwardMarksValidity(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<template<typename, typename> typename TSrCSA, typename TStorage,
    template<uint8_t> typename TAlphabet, uint8_t t_width, typename TBvValidMark>
void construct(SrCSAValidMark<TSrCSA<TStorage, TAlphabet<t_width>>, TBvValidMark> &t_index,
               const std::string &t_data_path,
               sri::Config &t_config) {
//  construct(dynamic_cast<TSrCSA &>(t_index), t_config);
  auto base_index = TSrCSA<TStorage, TAlphabet<t_width>>(t_index);
  construct(base_index, t_data_path, t_config);

  auto subsample_rate = t_index.SubsampleRate();
  auto prefix_key = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling validity marks and areas
    auto event = sdsl::memory_monitor::event("Subsampling Validity");
    auto key = prefix_key + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK;
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSubsamplingBackwardMarksValidity(subsample_rate, t_config);
    }

    if (!sdsl::cache_file_exists<TBvValidMark>(key, t_config)) {
      std::size_t r_prime;
      {
        sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(prefix_key + conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config));
        r_prime = buf.size();
      }
      constructBitVectorFromIntVector<TBvValidMark,
                                      typename TBvValidMark::rank_0_type,
                                      typename TBvValidMark::select_0_type>(key, t_config, r_prime, true);
    }
  }

  t_index.load(t_config);
}

template<typename TSrCSA, typename TBvValidMark>
void construct(SrCSAValidArea<TSrCSA, TBvValidMark> &t_index,
               const std::string &t_data_path,
               sri::Config &t_config) {
  construct(dynamic_cast<SrCSAValidMark<TSrCSA, TBvValidMark> &>(t_index), t_data_path, t_config);

  t_index.load(t_config);
}

inline void constructSubsamplingBackwardSamplesForPhiForwardWithBWTRuns(std::size_t t_subsample_rate,
                                                                        sdsl::cache_config& t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  // Samples
  sdsl::int_vector<> samples; // BWT-run starts positions in text
  sdsl::load_from_cache(samples, conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);
  auto r = samples.size();
  auto log_r = sdsl::bits::hi(r) + 1;

  sdsl::int_vector<> sorted_samples_idx;
  sdsl::load_from_cache(sorted_samples_idx, conf::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX, t_config);

  // We must sub-sample the samples associated to the first and last marks in the text
  const auto req_samples_idx = getExtremes(t_config, conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX);

  // Compute sub-sampled indices for sampled values
  sdsl::int_vector<> subsamples_idx;
  {
    auto subsamples_idx_vec = computeSampling(t_subsample_rate,
                                              std::reverse_iterator(sorted_samples_idx.end()),
                                              std::reverse_iterator(sorted_samples_idx.begin()),
                                              samples,
                                              req_samples_idx);

    std::sort(subsamples_idx_vec.begin(), subsamples_idx_vec.end());

    subsamples_idx = sdsl::int_vector<>(subsamples_idx_vec.size(), 0, log_r);
    std::copy(subsamples_idx_vec.begin(), subsamples_idx_vec.end(), subsamples_idx.begin());

    // Store sub-sample indices sorted by BWT positions
    sdsl::store_to_cache(subsamples_idx, prefix + conf::KEY_BWT_RUN_FIRST_IDX, t_config);
  }

  {
    // Compute sub-samples
    sdsl::int_vector<> subsamples(subsamples_idx.size(), 0, samples.width());
    std::transform(subsamples_idx.begin(), subsamples_idx.end(), subsamples.begin(),
                   [&samples](auto tt_i) { return samples[tt_i]; });

    sdsl::store_to_cache(subsamples, prefix + conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);
  }
}

template<uint8_t t_width>
auto computeSampleToMarkLinksForPhiForward(const std::string &t_prefix, sdsl::cache_config &t_config);

template<uint8_t t_width>
void constructSubsamplingBackwardMarksForPhiForwardWithBWTRuns(std::size_t t_subsample_rate,
                                                               sdsl::cache_config& t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  // Text positions of marks indices associated to sub-samples, i.e., text positions of sub-sampled marks
  sdsl::int_vector<> subsampled_mark_text_pos;
  std::size_t r_prime;

  {
    // Compute sub-sample to mark links

    // Links from sub-samples to their corresponding mark (actually, to rank of mark in BWT)
    auto subsample_to_mark_links = computeSampleToMarkLinksForPhiForward<t_width>(prefix, t_config);
    r_prime = subsample_to_mark_links.size();

    sdsl::int_vector<> bwt_run_ends_text_pos;
    sdsl::load_from_cache(bwt_run_ends_text_pos, conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config);

    subsampled_mark_text_pos = sdsl::int_vector(r_prime, 0, bwt_run_ends_text_pos.width());
    std::transform(subsample_to_mark_links.begin(),
                   subsample_to_mark_links.end(),
                   subsampled_mark_text_pos.begin(),
                   [&bwt_run_ends_text_pos](auto tt_i) { return bwt_run_ends_text_pos[tt_i]; });
  }

  auto log_r_prime = sdsl::bits::hi(r_prime) + 1;

  // Links from sub-sampled marks (sorted by text position) to sub-samples indices. Note that, initially, these are the indices of sub-sample in BWT.
  sdsl::int_vector<> subsampled_mark_to_subsample_links(r_prime, 0, log_r_prime);
  std::iota(subsampled_mark_to_subsample_links.begin(), subsampled_mark_to_subsample_links.end(), 0);

  // Sort indexes by text positions of its marks, becoming in the links from the sub-sampled marks to sub-sampled samples.
  std::sort(subsampled_mark_to_subsample_links.begin(),
            subsampled_mark_to_subsample_links.end(),
            [&subsampled_mark_text_pos](const auto &tt_a, const auto &tt_b) {
              return subsampled_mark_text_pos[tt_a] < subsampled_mark_text_pos[tt_b];
            });

  sdsl::store_to_cache(subsampled_mark_text_pos, prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config);

  sdsl::store_to_cache(subsampled_mark_to_subsample_links,
                       prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
                       t_config);
}

template<uint8_t t_width>
auto computeSampleToMarkLinksForPhiForward(const std::string &t_prefix, sdsl::cache_config &t_config) {
  // Sub-sampled indices of samples
  sdsl::int_vector<> subsamples_idx;
  sdsl::load_from_cache(subsamples_idx, t_prefix + conf::KEY_BWT_RUN_FIRST_IDX, t_config);

  sdsl::int_vector<> subsample_to_mark_links(subsamples_idx.size(), 0, subsamples_idx.width());

  // LF
  RLEString<> bwt_rle;
  sdsl::load_from_cache(bwt_rle, conf::KEY_BWT_RLE, t_config);
  auto get_char = buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = buildRankOfChar(std::cref(bwt_rle));

  Alphabet<t_width> alphabet;
  sdsl::load_from_cache(alphabet, conf::KEY_ALPHABET, t_config);
  auto n = alphabet.C[alphabet.sigma];

  auto get_f = [&alphabet](auto tt_symbol) { return alphabet.C[tt_symbol]; };
  auto lf = buildBasicLF(get_char, get_rank_of_char, get_f);

  // Psi
  PsiCoreRLE<> psi_core;
  sdsl::load_from_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };

  auto get_c = [&alphabet](auto tt_index) { return computeCForSAIndex(alphabet.C, tt_index); };
  auto cumulative = RandomAccessForCRefContainer(std::cref(alphabet.C));

  auto psi = Psi(psi_select, get_c, cumulative);

  // Marks positions
  sdsl::int_vector<> bwt_run_ends;
  sdsl::load_from_cache(bwt_run_ends, conf::KEY_BWT_RUN_LAST, t_config);

  auto rank_bwt_run_ends = [&bwt_run_ends](const auto &tt_k) {
    return std::lower_bound(bwt_run_ends.begin(), bwt_run_ends.end(), tt_k) - bwt_run_ends.begin();
  };

  // Samples positions
  sdsl::int_vector<> bwt_run_starts;
  sdsl::load_from_cache(bwt_run_starts, conf::KEY_BWT_RUN_FIRST, t_config);

  // Compute links from samples to marks
  for (int i = 0; i < subsamples_idx.size(); ++i) {
    subsample_to_mark_links[i] =
        computeSampleToMarkLinkForPhiForward(bwt_run_starts[subsamples_idx[i]], n, lf, psi, rank_bwt_run_ends);
  }

  return subsample_to_mark_links;
}

template<uint8_t t_width, typename TBvMark>
void constructSrCSACommonsWithBWTRuns(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  const auto n = sdsl::int_vector_buffer<>(cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config)).size();

  auto prefix = std::to_string(t_subsample_rate) + "_";

  // Sort samples (BWT-run last letter) by its text positions
  if (!sdsl::cache_file_exists(conf::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX, t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSortedIndices(conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config, conf::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX);
  }

  // Construct subsampling backward of samples (text positions of BWT-run first letter)
  if (!sdsl::cache_file_exists(prefix + conf::KEY_BWT_RUN_FIRST_IDX, t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardSamplesForPhiForwardWithBWTRuns(t_subsample_rate, t_config);
  }

  // Construct subsampling backward of marks (text positions of BWT-run last letter)
  if (!sdsl::cache_file_exists(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardMarksForPhiForwardWithBWTRuns<t_width>(t_subsample_rate, t_config);
  }

  // Construct successor on the text positions of sub-sampled BWT-run last letter
  if (!sdsl::cache_file_exists<TBvMark>(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config)) {
    auto event = sdsl::memory_monitor::event("Successor");
    constructBitVectorFromIntVector<TBvMark>(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config, n, false);
  }
}

inline void constructSamplesSortedByAlphabet(sdsl::cache_config &t_config) {
  PsiCoreRLE<> psi_rle;
  sdsl::load_from_cache(psi_rle, sdsl::conf::KEY_PSI, t_config, true);

  // Samples positions
  sdsl::int_vector<> samples;
  sdsl::load_from_cache(samples, conf::KEY_BWT_RUN_FIRST, t_config);

  const std::size_t buffer_size = 1 << 20;
  auto r = samples.size();
  auto log_r = sdsl::bits::hi(r) + 1;
  // BWT-run samples sorted by alphabet
  auto key = sdsl::cache_file_name(KeySortedByAlphabet(conf::KEY_BWT_RUN_FIRST_IDX), t_config);
  sdsl::int_vector_buffer<> samples_idx_sorted(key, std::ios::out, buffer_size, log_r);
  std::size_t n_runs = 0;
  auto report = [&samples, &samples_idx_sorted, &n_runs](const auto &tt_run_start, const auto &tt_run_end) {
    auto lower = std::lower_bound(samples.begin(), samples.end(), tt_run_start);
    samples_idx_sorted.push_back(lower - samples.begin());
    ++n_runs;
  };

  // Cumulative number of BWT-runs per symbols
  key = sdsl::cache_file_name(conf::KEY_BWT_RUN_CUMULATIVE_COUNT, t_config);
  sdsl::int_vector_buffer<> cumulative(key, std::ios::out, buffer_size, log_r);

  auto sigma = psi_rle.sigma();
  for (decltype(sigma) c = 0; c < sigma; ++c) {
    psi_rle.traverse(c, report);
    cumulative.push_back(n_runs);
  }

  samples_idx_sorted.close();
  sdsl::register_cache_file(KeySortedByAlphabet(conf::KEY_BWT_RUN_FIRST_IDX), t_config);
  cumulative.close();
  sdsl::register_cache_file(conf::KEY_BWT_RUN_CUMULATIVE_COUNT, t_config);
}

inline void constructSubsamplingBackwardSamplesSortedByAlphabet(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto key_prefix = std::to_string(t_subsample_rate) + "_";
  const std::size_t buffer_size = 1 << 20;

  sdsl::int_vector<> subsamples_idx_to_sorted;
  {
    sdsl::int_vector<> samples_idx_sorted;
    sdsl::load_from_cache(samples_idx_sorted, KeySortedByAlphabet(conf::KEY_BWT_RUN_FIRST_IDX), t_config);

    sdsl::bit_vector subsamples_idx_bv;
    {
      sdsl::int_vector<> subsamples_idx;
      sdsl::load_from_cache(subsamples_idx, key_prefix + conf::KEY_BWT_RUN_FIRST_IDX, t_config);

      subsamples_idx_bv = constructBitVectorFromIntVector(subsamples_idx, samples_idx_sorted.size(), false);

      subsamples_idx_to_sorted = sdsl::int_vector<>(subsamples_idx.size(), 0, subsamples_idx.width());
    }
    sdsl::bit_vector::rank_1_type rank_subsamples_idx_bv(&subsamples_idx_bv);

    sdsl::int_vector<> subsamples;
    sdsl::load_from_cache(subsamples, key_prefix + conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);

    auto key = sdsl::cache_file_name(KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_FIRST_TEXT_POS), t_config);
    sdsl::int_vector_buffer<> subsamples_sorted(key, std::ios::out, buffer_size, subsamples.width());

    key = sdsl::cache_file_name(KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_FIRST_IDX), t_config);
    sdsl::int_vector_buffer<> subsamples_idx_sorted(key, std::ios::out, buffer_size, subsamples_idx_to_sorted.width());

    std::size_t i = 0;
    std::size_t c_subsamples = 0;
    for (auto &&idx : samples_idx_sorted) {
      if (subsamples_idx_bv[idx] == true) {
        auto rnk = rank_subsamples_idx_bv(idx);
        subsamples_sorted.push_back(subsamples[rnk]);
        subsamples_idx_sorted.push_back(i);

        subsamples_idx_to_sorted[rnk] = c_subsamples++;
      }
      ++i;
    }

    subsamples_sorted.close();
    sdsl::register_cache_file(KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_FIRST_TEXT_POS), t_config);
    subsamples_idx_sorted.close();
    sdsl::register_cache_file(KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_FIRST_IDX), t_config);
  }

  sdsl::int_vector<> mark_to_sample_idx;
  sdsl::load_from_cache(mark_to_sample_idx, key_prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config);
  auto key = sdsl::cache_file_name(
      KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX), t_config);
  sdsl::int_vector_buffer<> mark_to_sample_idx_sorted(key, std::ios::out, buffer_size, mark_to_sample_idx.width());
  for (const auto &idx : mark_to_sample_idx) {
    mark_to_sample_idx_sorted.push_back(subsamples_idx_to_sorted[idx]);
  }

  mark_to_sample_idx_sorted.close();
  sdsl::register_cache_file(KeySortedByAlphabet(key_prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX),
                            t_config);
}

inline void constructSubsamplingBackwardSamplesPosition(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  sdsl::int_vector<> samples_pos; // BWT-run starts positions in SA
  sdsl::load_from_cache(samples_pos, conf::KEY_BWT_RUN_FIRST, t_config);

  sdsl::int_vector<> subsamples_idx; // Sub-samples indices
  sdsl::load_from_cache(subsamples_idx, prefix + conf::KEY_BWT_RUN_FIRST_IDX, t_config);

  // Compute sub-samples positions
  sdsl::int_vector<> subsamples_pos(subsamples_idx.size(), 0, samples_pos.width());
  std::transform(subsamples_idx.begin(), subsamples_idx.end(), subsamples_pos.begin(),
                 [&samples_pos](auto tt_i) { return samples_pos[tt_i]; });

  sdsl::store_to_cache(subsamples_pos, prefix + conf::KEY_BWT_RUN_FIRST, t_config);
}

inline void constructSubsamplingBackwardMarksValidity(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  sdsl::int_vector<> marks;
  sdsl::load_from_cache(marks, conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config);
  sdsl::int_vector<> sorted_marks_idx;
  sdsl::load_from_cache(sorted_marks_idx, conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX, t_config);
  auto it_marks_idx = sorted_marks_idx.end();
  auto get_next_mark = [&marks, &it_marks_idx]() { return marks[*(--it_marks_idx)]; };

  sdsl::int_vector<> submarks;
  sdsl::load_from_cache(submarks, prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config);
  std::sort(submarks.begin(), submarks.end());
  auto it_submarks = submarks.end();
  auto get_next_submark = [&it_submarks]() { return *(--it_submarks); };

  auto r_prime = submarks.size();
  std::vector<std::pair<std::size_t, std::size_t>> validity;
  validity.reserve(r_prime / 4);
  std::size_t max_valid_area = 0;
  auto report = [r_prime, &validity, &max_valid_area](auto tt_i, auto tt_submark, auto tt_next_mark) {
    auto valid_area = tt_submark - tt_next_mark;
    validity.emplace_back(r_prime - tt_i - 1, valid_area);

    if (max_valid_area < valid_area) max_valid_area = valid_area;
  };

  computeSubmarksValidity(r_prime, get_next_mark, get_next_submark, report);

  const std::size_t buffer_size = 1 << 20;
  sdsl::int_vector_buffer<> valid_submarks(
      sdsl::cache_file_name(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK, t_config),
      std::ios::out,
      buffer_size,
      sdsl::bits::hi(r_prime) + 1);

  sdsl::int_vector_buffer<> valid_areas(
      sdsl::cache_file_name(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA, t_config),
      std::ios::out,
      buffer_size,
      sdsl::bits::hi(max_valid_area) + 1);

  for (auto it = validity.rbegin(); it != validity.rend(); ++it) {
    valid_submarks.push_back(it->first);
    valid_areas.push_back(it->second);
  }

  valid_submarks.close();
  sdsl::register_cache_file(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK, t_config);
  valid_areas.close();
  sdsl::register_cache_file(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA, t_config);
}

}

#endif //SRI_SR_CSA_H_
