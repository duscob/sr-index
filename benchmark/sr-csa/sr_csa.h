//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_SR_CSA_H_
#define SRI_BENCHMARK_SR_CSA_SR_CSA_H_

#include <optional>

#include "sr-index/sampling.h"

#include "csa.h"

struct RunCSA {
  std::size_t start = 0;
  std::size_t end = 0;
};

auto limits(const RunCSA &t_run) {
  return t_run;
}

auto &limits(RunCSA &t_run) {
  return t_run;
}

bool operator<(const RunCSA &lhs, const RunCSA &rhs) {
  return lhs.start < rhs.start;
}

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSampleIdx = sdsl::sd_vector<>,
    typename TRun = RunCSA>
class SrCSA : public CSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample> {
 public:
  using BaseClass = CSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample>;

  explicit SrCSA(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr)
      : BaseClass(t_storage), subsample_rate_{t_sr}, key_prefix_{std::to_string(subsample_rate_) + "_"} {
  }

  virtual std::size_t SubsampleRate() const { return subsample_rate_; }

  using typename BaseClass::size_type;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(key(SrIndexKey::ALPHABET), out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), out, child, "psi");

    written_bytes += this->template serializeItem<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_select");

    written_bytes +=
        this->template serializeItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), out, child, "mark_to_sample");

    written_bytes += this->template serializeItem<TSample>(key(SrIndexKey::SAMPLES), out, child, "samples");

    written_bytes +=
        this->template serializeItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), out, child, "samples_idx");

    return written_bytes;
  }

 protected:

  using BaseClass::key;
  void setupKeyNames() override {
    if (!this->keys_.empty()) return;

    BaseClass::setupKeyNames();
    this->keys_.resize(6);
//    key(SrIndexKey::ALPHABET) = sri::key_trait<t_width>::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
    key(SrIndexKey::MARKS) = key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
    key(SrIndexKey::SAMPLES) = key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS;
    key(SrIndexKey::MARK_TO_SAMPLE) =
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
    key(SrIndexKey::SAMPLES_IDX) = key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX;
  }

  using typename BaseClass::TSource;

  using typename BaseClass::Value;
  using TFnGetSample = std::function<std::optional<Value>(Value)>;
  virtual TFnGetSample constructGetSampleForSAIdx(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto get_run = [cref_psi_core](auto tt_sa_pos) {
      auto[run, run_start] = cref_psi_core.get().rankSoftRun(tt_sa_pos);
      return std::make_pair(run - 1, run_start);
    };

    auto cref_bv_sample_idx = this->template loadItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto is_run_sampled = [cref_bv_sample_idx](auto tt_i) {
      return cref_bv_sample_idx.get()[tt_i];
    };

    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto bv_sample_idx_rank = this->template loadBVRank<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto get_sample = [cref_samples, bv_sample_idx_rank](auto tt_run) {
      auto sampled_run = bv_sample_idx_rank(tt_run);
      return cref_samples.get()[sampled_run];
    };

    return sri::GetSampleForSAPosition(get_run, is_run_sampled, get_sample);
  }

  using typename BaseClass::RunData;
  using TFnGetSampleForRunData = std::function<std::optional<Value>(const RunData *)>;
  virtual TFnGetSampleForRunData constructGetSampleForRunData(TSource &t_source) {
    auto get_sample_for_sa_idx = constructGetSampleForSAIdx(t_source);

    return [get_sample_for_sa_idx](const RunData *tt_run_data) {
      return get_sample_for_sa_idx(tt_run_data->pos);
    };
  }

  using TFnGetSampleForRun = std::function<std::optional<Value>(const TRun &)>;
  virtual TFnGetSampleForRun constructGetSampleForRun(TSource &t_source) {
    auto get_sample_for_sa_idx = constructGetSampleForSAIdx(t_source);

    return [get_sample_for_sa_idx](const TRun &tt_run) {
      return get_sample_for_sa_idx(tt_run.start);
    };
  }

  using Char = typename TPsiRLE::Char;
  using TFnCreateRun = std::function<TRun(std::size_t, std::size_t, Char, std::size_t, bool)>;
  virtual TFnCreateRun constructCreateRun() {
    auto create_run = [](auto tt_first, auto tt_last, auto, auto, auto) {
      return TRun{tt_first, tt_last};
    };

    return create_run;
  }

  using typename BaseClass::Range;
  using Runs = std::vector<TRun>;
  using typename BaseClass::Position;
  using TFnSplitRange = std::function<Runs(Range)>;
  virtual TFnSplitRange constructSplitRangeInBWTRuns(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto create_run = constructCreateRun();

    auto split = [cref_psi_core, create_run](const auto &tt_range) -> auto {
      const auto &[first, last] = tt_range;
      return cref_psi_core.get().splitInSortedRuns(first, last, create_run);
    };

    return split;
  }

  using TFnSplitRun = std::function<Runs(TRun)>;
  virtual TFnSplitRun constructSplitRunInBWTRuns(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto cref_alphabet = this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);
    auto create_run = constructCreateRun();

    auto split = [cref_psi_core, cref_alphabet, create_run](const auto &tt_run) -> auto {
      const auto &[first, last] = limits(tt_run);
      const auto &cumulative = cref_alphabet.get().C;
      auto c = sri::computeCForSAIndex(cumulative, first);

      auto cum_c = cumulative[c];
      auto cum_next_c = cumulative[c + 1];
      auto first_rank = first - cum_c + 1;
      auto last_rank = std::min(last, cum_next_c) - cum_c + 1;
      Runs runs;
      auto report = [&runs, &create_run](auto tt_first, auto tt_last, auto tt_c, auto tt_n_run, auto tt_is_first) {
        runs.emplace_back(create_run(tt_first, tt_last, tt_c, tt_n_run, tt_is_first));
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

    return split;
  }

  using BvMark = TBvMark;
  using MarkToSampleIdx = TMarkToSampleIdx;
  using Sample = TSample;
  using Run = TRun;
  using typename BaseClass::TFnPhiForRange;
  TFnPhiForRange constructPhiForRange(TSource &t_source) override {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto successor = sri::CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), t_source);
    auto get_mark_to_sample_idx = sri::RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, false);

    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto get_sample = sri::RandomAccessForCRefContainer(cref_samples);
    sri::SampleValidatorDefault sample_validator_default;

    auto phi = sri::buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);
    auto phi_simple = [phi](const auto &tt_prev_value) { return phi(tt_prev_value).first; };

    auto get_sample_4_run = constructGetSampleForRun(t_source);
    auto split_range = constructSplitRangeInBWTRuns(t_source);
    auto split_run = constructSplitRunInBWTRuns(t_source);

    auto is_range_empty = this->constructIsRangeEmpty();

    auto update_run = [](TRun &tt_run) {
      auto &[first, last] = limits(tt_run);
      ++first;
      return tt_run;
    };

    auto is_run_empty = [](const TRun &tt_run) {
      const auto &[first, last] = limits(tt_run);
      return !(first < last);
    };

    return sri::PhiForwardForRange(phi_simple,
                                   get_sample_4_run,
                                   split_range,
                                   split_run,
                                   subsample_rate_,
                                   this->n_,
                                   is_range_empty,
                                   update_run,
                                   is_run_empty);
  }

  using TFnPsiForRunData = std::function<RunData * (RunData * )>;
  virtual TFnPsiForRunData constructPsiForRunData(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto psi_select = [cref_psi_core](auto tt_c, auto tt_rnk) { return cref_psi_core.get().select(tt_c, tt_rnk); };

    auto cref_alphabet = this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);
    auto get_c = [cref_alphabet](auto tt_index) { return sri::computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = sri::RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));

    auto psi = sri::Psi(psi_select, get_c, cumulative);

    return [psi](RunData *tt_run_data) {
      tt_run_data->pos = psi(tt_run_data->pos);
      return tt_run_data;
    };
  }

  using typename BaseClass::TFnComputeToehold;
  TFnComputeToehold constructComputeToehold(TSource &t_source) override {
    auto get_sample = constructGetSampleForRunData(t_source);
    auto psi_4_run_data = constructPsiForRunData(t_source);
    auto compute_sa_value_for_bwt_run_start = sri::buildComputeSAValueForward(get_sample, psi_4_run_data, this->n_);
    auto compute_sa_value = [compute_sa_value_for_bwt_run_start](const std::shared_ptr<RunData> &tt_run_data) {
      return compute_sa_value_for_bwt_run_start(tt_run_data.get());
    };

    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);

    return sri::buildComputeToeholdForPhiForward(compute_sa_value, cref_psi_core.get().size());
  }

  std::size_t subsample_rate_ = 1;
  std::string key_prefix_;
};

template<typename TChar>
struct RunSrCSA : RunCSA {
  TChar c;
  std::size_t partial_rank;
  bool is_run_start;
};

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSampleIdx = sdsl::sd_vector<>,
    typename TRunCumulativeCount = sdsl::int_vector<>,
    typename TRun = RunSrCSA<typename TPsiRLE::Char>>
class SrCSASlim : public SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSampleIdx, TRun> {
 public:
  using BaseClass = SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSampleIdx, TRun>;

  SrCSASlim(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr) : BaseClass(t_storage, t_sr) {}

  using typename BaseClass::size_type;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(key(SrIndexKey::ALPHABET), out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), out, child, "psi");

    written_bytes += this->template serializeItem<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(key(SrIndexKey::MARKS), out, child, "marks_select");

    written_bytes +=
        this->template serializeItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), out, child, "mark_to_sample");

    written_bytes += this->template serializeItem<TSample>(key(SrIndexKey::SAMPLES), out, child, "samples");

    written_bytes +=
        this->template serializeItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), out, child, "samples_idx");

    written_bytes += this->template serializeItem<TRunCumulativeCount>(
        key(SrIndexKey::RUN_CUMULATIVE_COUNT), out, child, "run_cumulative_count");

    return written_bytes;
  }

 protected:

  using BaseClass::key;
  void setupKeyNames() override {
    if (!this->keys_.empty()) return;

    BaseClass::setupKeyNames();
    this->keys_.resize(7);
//    key(SrIndexKey::ALPHABET) = sri::key_trait<t_width>::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
    key(SrIndexKey::MARKS) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
    key(SrIndexKey::SAMPLES) = sri::KeySortedByAlphabet(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS);
    key(SrIndexKey::MARK_TO_SAMPLE) = sri::KeySortedByAlphabet(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX);
    key(SrIndexKey::SAMPLES_IDX) = sri::KeySortedByAlphabet(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX);
    key(SrIndexKey::RUN_CUMULATIVE_COUNT) = sri::key_trait<t_width>::KEY_BWT_RUN_CUMULATIVE_COUNT;
  }

  using typename BaseClass::TSource;

  using typename BaseClass::RunData;
  using typename BaseClass::Char;
  struct RunDataExt : RunData {
    Char c;
    std::size_t partial_rank;
    bool is_run_start;

    RunDataExt(std::size_t t_pos = 0, Char t_c = 0, std::size_t t_partial_rank = 0, bool t_is_run_start = true)
        : RunData(t_pos), c{t_c}, partial_rank{t_partial_rank}, is_run_start{t_is_run_start} {}
  };

  using typename BaseClass::Range;
  using typename BaseClass::RangeLF;
  using typename BaseClass::DataBackwardSearchStep;
  using typename BaseClass::TFnCreateDataBackwardSearchStep;
  TFnCreateDataBackwardSearchStep constructCreateDataBackwardSearchStep() override {
    return [](const Range &tt_range, Char tt_c, const RangeLF &tt_next_range, std::size_t tt_step) {
      const auto &[start, end] = tt_next_range;
      return DataBackwardSearchStep{tt_step,
                                    std::make_shared<RunDataExt>(
                                        start.run.start, tt_c, start.run.rank, tt_range.start != start.run.start)};
    };
  }

  using typename BaseClass::TFnGetInitialDataBackwardSearchStep;
  TFnGetInitialDataBackwardSearchStep constructGetInitialDataBackwardSearchStep(TSource &t_source) override {
    return [](const auto &tt_step) { return DataBackwardSearchStep{0, std::make_shared<RunDataExt>()}; };
  }

  using typename BaseClass::Value;
  auto constructGetSample(TSource &t_source) {
    auto cref_run_cum_c = this->template loadItem<TRunCumulativeCount>(key(SrIndexKey::RUN_CUMULATIVE_COUNT), t_source);
    auto cref_bv_sample_idx = this->template loadItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_idx_rank = this->template loadBVRank<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);

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

  using typename BaseClass::TFnGetSampleForRunData;
  virtual TFnGetSampleForRunData constructGetSampleForRunData(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);

    return [get_sample](const RunData *tt_run_data) {
      auto run_data_ext = dynamic_cast<const RunDataExt *>(tt_run_data);
      return get_sample(run_data_ext->c, run_data_ext->partial_rank, run_data_ext->is_run_start);
    };
  }

  using typename BaseClass::TFnGetSampleForRun;
  virtual TFnGetSampleForRun constructGetSampleForRun(TSource &t_source) {
    auto get_sample = constructGetSample(t_source);

    return [get_sample](const TRun &tt_run) {
      return get_sample(tt_run.c, tt_run.partial_rank, tt_run.is_run_start);
    };
  }

  using typename BaseClass::TFnCreateRun;
  TFnCreateRun constructCreateRun() override {
    auto create_run = [](auto tt_first, auto tt_last, auto tt_c, auto tt_partial_rank, auto tt_is_first) {
      return TRun{tt_first, tt_last, tt_c, tt_partial_rank, tt_is_first};
    };

    return create_run;
  }

  using typename BaseClass::TFnPsiForRunData;
  TFnPsiForRunData constructPsiForRunData(TSource &t_source) override {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto psi_select = [cref_psi_core](Char tt_c, auto tt_rnk) {
      RunDataExt run_data;
      auto report = [&run_data, &tt_c](auto tt_pos, auto tt_run_start, auto tt_run_end, auto tt_run_rank) {
        run_data = RunDataExt{tt_pos, tt_c, tt_run_rank, tt_pos == tt_run_start};
      };
      cref_psi_core.get().select(tt_c, tt_rnk, report);
      return run_data;
    };

    auto cref_alphabet = this->template loadItem<TAlphabet>(key(SrIndexKey::ALPHABET), t_source);
    auto get_c = [cref_alphabet](auto tt_index) { return sri::computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = sri::RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));

    auto psi = sri::Psi(psi_select, get_c, cumulative);

    return [psi](RunData *tt_run_data) {
      auto run_data_ext = dynamic_cast<RunDataExt *>(tt_run_data);
      *run_data_ext = psi(tt_run_data->pos);
      return run_data_ext;
    };
  }
};

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSamplePos = sdsl::sd_vector<>,
    typename TRun = RunCSA>
class SrCSAWithBv : public SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos, TRun> {
 public:
  using BaseClass = SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos, TRun>;

  SrCSAWithBv(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr) : BaseClass(t_storage, t_sr) {}

 protected:

  using BaseClass::key;
  void setupKeyNames() override {
    if (!this->keys_.empty()) return;

    BaseClass::setupKeyNames();
//    this->keys_.resize(6);
//    key(SrIndexKey::ALPHABET) = sri::key_trait<t_width>::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
//    key(SrIndexKey::MARKS) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
//    key(SrIndexKey::SAMPLES) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS;
//    key(SrIndexKey::MARK_TO_SAMPLE) =
//        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
    key(SrIndexKey::SAMPLES_IDX) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST;
  }

  using typename BaseClass::TSource;

  using typename BaseClass::TFnGetSample;
  TFnGetSample constructGetSampleForSAIdx(TSource &t_source) override {
    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto cref_bv_sample_pos = this->template loadItem<TBvSamplePos>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(key(SrIndexKey::SAMPLES_IDX), t_source, true);

    auto get_sample = [cref_samples, cref_bv_sample_pos, bv_sample_pos_rank](auto tt_pos) {
      return cref_bv_sample_pos.get()[tt_pos]
             ? std::optional<std::size_t>(cref_samples.get()[bv_sample_pos_rank(tt_pos)])
             : std::nullopt;
    };

    return get_sample;
  }

  using typename BaseClass::Runs;
  using typename BaseClass::TFnSplitRange;
  TFnSplitRange constructSplitRangeInBWTRuns(TSource &t_source) override {
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_pos_select = this->template loadBVSelect<TBvSamplePos>(key(SrIndexKey::SAMPLES_IDX), t_source, true);

    auto split = [bv_sample_pos_rank, bv_sample_pos_select](const auto &tt_range) -> auto {
      const auto &[first, last] = tt_range;
      auto first_rank = bv_sample_pos_rank(first + 1) + 1;
      auto last_rank = bv_sample_pos_rank(last);

      Runs runs;
      auto run_start = first;
      for (int i = first_rank; i <= last_rank; ++i) {
        auto run_end = bv_sample_pos_select(i);
        runs.emplace_back(TRun{run_start, run_end});

        run_start = run_end;
      }

      runs.emplace_back(TRun{run_start, last});

      return runs;
    };

    return split;
  }
};

template<typename TSrCSA, typename TBvValidMark = sdsl::bit_vector>
class SrCSAValidMark : public TSrCSA {
 public:
  using BaseClass = TSrCSA;

  SrCSAValidMark(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr) : BaseClass(t_storage, t_sr) {}

  using typename BaseClass::size_type;
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = BaseClass::serialize(out, v, name);
    written_bytes +=
        this->template serializeItem<TBvValidMark>(key(SrIndexKey::VALID_MARKS), out, child, "valid_marks");

    return written_bytes;
  }

 protected:

  using BaseClass::key;
  void setupKeyNames() override {
    if (!this->keys_.empty()) return;

    BaseClass::setupKeyNames();
    this->keys_.resize(8);
    constexpr uint8_t t_width = SrCSAValidMark<TSrCSA>::AlphabetWidth;
//    key(SrIndexKey::ALPHABET) = sri::key_trait<t_width>::KEY_ALPHABET;
//    key(SrIndexKey::NAVIGATE) = sdsl::conf::KEY_PSI;
//    key(SrIndexKey::MARKS) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
//    key(SrIndexKey::SAMPLES) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS;
//    key(SrIndexKey::MARK_TO_SAMPLE) =
//        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;
//    key(SrIndexKey::SAMPLES_IDX) = this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX;
    key(SrIndexKey::VALID_MARKS) =
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK;
  }

  using typename BaseClass::TSource;

  using typename BaseClass::BvMark;
  using typename BaseClass::MarkToSampleIdx;
  using typename BaseClass::Sample;
  using typename BaseClass::Run;
  using typename BaseClass::TFnPhiForRange;
  TFnPhiForRange constructPhiForRange(TSource &t_source) override {
    auto bv_mark_rank = this->template loadBVRank<BvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<BvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto successor = sri::CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<MarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), t_source);
    auto cref_valid_mark = this->template loadItem<TBvValidMark>(key(SrIndexKey::VALID_MARKS), t_source, true);
    auto get_mark_to_sample_idx = sri::RandomAccessForTwoContainers(cref_mark_to_sample_idx, cref_valid_mark);

    auto cref_samples = this->template loadItem<Sample>(key(SrIndexKey::SAMPLES), t_source);
    auto get_sample = sri::RandomAccessForCRefContainer(cref_samples);
    sri::SampleValidatorDefault sample_validator_default;

    auto phi = sri::buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);

    auto get_sample_4_run = this->constructGetSampleForRun(t_source);
    auto split_range = this->constructSplitRangeInBWTRuns(t_source);
    auto split_run = this->constructSplitRunInBWTRuns(t_source);

    auto is_range_empty = this->constructIsRangeEmpty();

    auto update_run = [](Run &tt_run) {
      auto &[first, last] = limits(tt_run);
      ++first;
      return tt_run;
    };

    auto is_run_empty = [](const Run &tt_run) {
      const auto &[first, last] = limits(tt_run);
      return !(first < last);
    };

    return sri::PhiForwardForRangeWithValidity(phi,
                                               get_sample_4_run,
                                               split_range,
                                               split_run,
                                               this->subsample_rate_,
                                               this->n_,
                                               is_range_empty,
                                               update_run,
                                               is_run_empty);
  }
};

template<uint8_t t_width, typename TBVMark>
void constructSrCSACommons(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBVMark, typename TMarkToSampleIdx, typename TSample, typename TBVSampleIdx, typename TRun>
void construct(SrCSA<t_width, TAlphabet, TPsiCore, TBVMark, TMarkToSampleIdx, TSample, TBVSampleIdx, TRun> &t_index,
               sdsl::cache_config &t_config) {
  auto subsample_rate = t_index.SubsampleRate();

  constructSrCSACommons<t_width, TBVMark>(subsample_rate, t_config);

  std::size_t r;
  {
    sdsl::int_vector_buffer<> bwt(sdsl::cache_file_name(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config));
    r = bwt.size();
  }

  auto prefix = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling backward of samples (text positions of BWT-run last letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX;
    if (!sdsl::cache_file_exists<TBVSampleIdx>(key, t_config)) {
      sri::constructBitVectorFromIntVector<TBVSampleIdx>(key, t_config, r, false);
    }
  }

  t_index.load(t_config);
}

template<uint8_t t_width>
void constructSamplesSortedByAlphabet(sdsl::cache_config &t_config);

template<uint8_t t_width>
void constructSubsamplingBackwardSamplesSortedByAlphabet(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBVMark, typename TMarkToSample, typename TSample, typename TBVSampleIdx, typename TRunCumCnt, typename TRun>
void construct(
    SrCSASlim<t_width, TAlphabet, TPsiCore, TBVMark, TMarkToSample, TSample, TBVSampleIdx, TRunCumCnt, TRun> &t_index,
    sdsl::cache_config &t_config) {
  auto subsample_rate = t_index.SubsampleRate();

  constructSrCSACommons<t_width, TBVMark>(subsample_rate, t_config);

  {
    // Construct samples' indices sorted by alphabet
    auto event = sdsl::memory_monitor::event("Samples");
    auto key = sri::KeySortedByAlphabet(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX);
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSamplesSortedByAlphabet<t_width>(t_config);
    }
  }

  auto prefix_key = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling backward of samples sorted by alphabet
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = sri::KeySortedByAlphabet(prefix_key + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS);
    if (!sdsl::cache_file_exists<TBVSampleIdx>(key, t_config)) {
      constructSubsamplingBackwardSamplesSortedByAlphabet<t_width>(subsample_rate, t_config);
    }
  }

  {
    // Construct subsampling indices backward of samples sorted by alphabet
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = sri::KeySortedByAlphabet(prefix_key + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX);
    if (!sdsl::cache_file_exists<TBVSampleIdx>(key, t_config)) {
      std::size_t r;
      {
        sdsl::int_vector_buffer<> bwt(sdsl::cache_file_name(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config));
        r = bwt.size();
      }

      sri::constructBitVectorFromIntVector<TBVSampleIdx>(key, t_config, r, false);
    }
  }

  t_index.load(t_config);
}

template<uint8_t t_width>
void constructSubsamplingBackwardSamplesPosition(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBvMark, typename TMarkToSampleIdx, typename TSample, typename TBvSamplePos>
void construct(SrCSAWithBv<t_width, TAlphabet, TPsiCore, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos> &t_index,
               sdsl::cache_config &t_config) {
  auto subsample_rate = t_index.SubsampleRate();

  constructSrCSACommons<t_width, TBvMark>(subsample_rate, t_config);

  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = buf.size();
  }

  auto prefix = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling backward of samples (text positions of BWT-run last letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST;
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSubsamplingBackwardSamplesPosition<t_width>(subsample_rate, t_config);
    }

    if (!sdsl::cache_file_exists<TBvSamplePos>(key, t_config)) {
      sri::constructBitVectorFromIntVector<TBvSamplePos>(key, t_config, n, false);
    }
  }

  t_index.load(t_config);
}

template<uint8_t t_width>
void constructSubsamplingBackwardMarksValidity(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<typename TSrCSA, typename TBvValidMark>
void construct(SrCSAValidMark<TSrCSA, TBvValidMark> &t_index, sdsl::cache_config &t_config) {
//  construct(dynamic_cast<TSrCSA &>(t_index), t_config);
  auto base_index = TSrCSA(t_index);
  construct(base_index, t_config);

  auto subsample_rate = t_index.SubsampleRate();
  auto prefix_key = std::to_string(subsample_rate) + "_";

  {
    // Construct subsampling validity marks and areas
    auto event = sdsl::memory_monitor::event("Subsampling Validity");
    constexpr uint8_t t_width = SrCSAValidMark<TSrCSA>::AlphabetWidth;
    auto key = prefix_key + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK;
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSubsamplingBackwardMarksValidity<t_width>(subsample_rate, t_config);
    }

    if (!sdsl::cache_file_exists<TBvValidMark>(key, t_config)) {
      std::size_t r_prime;
      {
        sdsl::int_vector_buffer<> buf(sdsl::cache_file_name(
            prefix_key + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config));
        r_prime = buf.size();
      }
      sri::constructBitVectorFromIntVector<TBvValidMark>(key, t_config, r_prime, true);
    }
  }

  t_index.load(t_config);
}

template<uint8_t t_width>
void constructSubsamplingBackwardSamples(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  // Samples
  sdsl::int_vector<> samples; // BWT-run starts positions in text
  sdsl::load_from_cache(samples, sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);
  auto r = samples.size();
  auto log_r = sdsl::bits::hi(r) + 1;

  sdsl::int_vector<> sorted_samples_idx;
  sdsl::load_from_cache(sorted_samples_idx, sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX, t_config);

  std::array<std::size_t, 2> req_samples_idx{};
  {
    // We must sub-sample the samples associated to the first and last marks in the text
    sdsl::int_vector_buffer<> mark_to_sample(
        sdsl::cache_file_name(sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config));
    req_samples_idx[0] = mark_to_sample[0];
    req_samples_idx[1] = mark_to_sample[mark_to_sample.size() - 1];
  }

  // Compute sub-sampled indices for sampled values
  sdsl::int_vector<> subsamples_idx;
  {
    auto subsamples_idx_vec = sri::computeSampling(t_subsample_rate,
                                                   std::reverse_iterator(sorted_samples_idx.end()),
                                                   std::reverse_iterator(sorted_samples_idx.begin()),
                                                   samples,
                                                   req_samples_idx);

    std::sort(subsamples_idx_vec.begin(), subsamples_idx_vec.end());

    subsamples_idx = sdsl::int_vector<>(subsamples_idx_vec.size(), 0, log_r);
    std::copy(subsamples_idx_vec.begin(), subsamples_idx_vec.end(), subsamples_idx.begin());

    // Store sub-sample indices sorted by BWT positions
    sdsl::store_to_cache(subsamples_idx, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_config);
  }

  {
    // Compute sub-samples
    sdsl::int_vector<> subsamples(subsamples_idx.size(), 0, samples.width());
    std::transform(subsamples_idx.begin(), subsamples_idx.end(), subsamples.begin(),
                   [&samples](auto tt_i) { return samples[tt_i]; });

    sdsl::store_to_cache(subsamples, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);
  }
}

template<uint8_t t_width>
auto computeSampleToMarkLinks(const std::string &t_prefix, sdsl::cache_config &t_config);

template<uint8_t t_width>
void constructSubsamplingBackwardMarks(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  // Text positions of marks indices associated to sub-samples, i.e., text positions of sub-sampled marks
  sdsl::int_vector<> subsampled_mark_text_pos;
  std::size_t r_prime;

  {
    // Compute sub-sample to mark links

    // Links from sub-samples to their corresponding mark (actually, to rank of mark in BWT)
    auto subsample_to_mark_links = computeSampleToMarkLinks<t_width>(prefix, t_config);
    r_prime = subsample_to_mark_links.size();

    sdsl::int_vector<> bwt_run_ends_text_pos;
    sdsl::load_from_cache(bwt_run_ends_text_pos, sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_config);

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

  sdsl::store_to_cache(subsampled_mark_text_pos,
                       prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST,
                       t_config);

  sdsl::store_to_cache(subsampled_mark_to_subsample_links,
                       prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
                       t_config);
}

template<uint8_t t_width>
auto computeSampleToMarkLinks(const std::string &t_prefix, sdsl::cache_config &t_config) {
  // Sub-sampled indices of samples
  sdsl::int_vector<> subsamples_idx;
  sdsl::load_from_cache(subsamples_idx, t_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_config);

  sdsl::int_vector<> subsample_to_mark_links(subsamples_idx.size(), 0, subsamples_idx.width());

  // LF
  sri::rle_string<> bwt_rle;
  sdsl::load_from_cache(bwt_rle, sri::key_trait<t_width>::KEY_BWT_RLE, t_config);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  typename sri::alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, sri::key_trait<t_width>::KEY_ALPHABET, t_config);
  auto n = alphabet.C[alphabet.sigma];

  auto get_f = [&alphabet](auto tt_symbol) { return alphabet.C[alphabet.char2comp[tt_symbol]]; };
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  // Psi
  sri::PsiCoreRLE<> psi_core;
  sdsl::load_from_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };

  auto get_c = [&alphabet](auto tt_index) { return sri::computeCForSAIndex(alphabet.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  // Marks positions
  sdsl::int_vector<> bwt_run_ends;
  sdsl::load_from_cache(bwt_run_ends, sri::key_trait<t_width>::KEY_BWT_RUN_LAST, t_config);

  auto rank_bwt_run_ends = [&bwt_run_ends](const auto &tt_k) {
    return std::lower_bound(bwt_run_ends.begin(), bwt_run_ends.end(), tt_k) - bwt_run_ends.begin();
  };

  // Samples positions
  sdsl::int_vector<> bwt_run_starts;
  sdsl::load_from_cache(bwt_run_starts, sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);

  // Compute links from samples to marks
  for (int i = 0; i < subsamples_idx.size(); ++i) {
    subsample_to_mark_links[i] =
        sri::computeSampleToMarkLinkForPhiForward(bwt_run_starts[subsamples_idx[i]], n, lf, psi, rank_bwt_run_ends);
  }

  return subsample_to_mark_links;
}

template<uint8_t t_width, typename TBvMark>
void constructSrCSACommons(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = bwt_buf.size();
  }

  auto prefix = std::to_string(t_subsample_rate) + "_";

  {
    // Sort samples (BWT-run last letter) by its text positions
    auto event = sdsl::memory_monitor::event("Subsampling");
    const auto key = sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX;
    if (!sdsl::cache_file_exists<TBvMark>(key, t_config)) {
      sri::constructSortedIndices(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config, key);
    }
  }

  {
    // Construct subsampling backward of samples (text positions of BWT-run last letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    auto key = prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX;
    if (!sdsl::cache_file_exists(key, t_config)) {
      constructSubsamplingBackwardSamples<t_width>(t_subsample_rate, t_config);
    }
  }

  {
    // Construct subsampling backward of marks (text positions of BWT-run first letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    if (!sdsl::cache_file_exists(
        prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config)) {
      constructSubsamplingBackwardMarks<t_width>(t_subsample_rate, t_config);
    }
  }

  {
    // Construct successor on the text positions of sub-sampled BWT-run last letter
    auto event = sdsl::memory_monitor::event("Successor");
    const auto key = prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST;
    if (!sdsl::cache_file_exists<TBvMark>(key, t_config)) {
      sri::constructBitVectorFromIntVector<TBvMark>(key, t_config, n, false);
    }
  }
}

template<uint8_t t_width>
void constructSamplesSortedByAlphabet(sdsl::cache_config &t_config) {
  sri::PsiCoreRLE<> psi_rle;
  sdsl::load_from_cache(psi_rle, sdsl::conf::KEY_PSI, t_config, true);

  // Samples positions
  sdsl::int_vector<> samples;
  sdsl::load_from_cache(samples, sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);

  const std::size_t buffer_size = 1 << 20;
  auto r = samples.size();
  auto log_r = sdsl::bits::hi(r) + 1;
  // BWT-run samples sorted by alphabet
  auto key = sdsl::cache_file_name(sri::KeySortedByAlphabet(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX), t_config);
  sdsl::int_vector_buffer<> samples_idx_sorted(key, std::ios::out, buffer_size, log_r);
  std::size_t n_runs = 0;
  auto report = [&samples, &samples_idx_sorted, &n_runs](const auto &tt_run_start, const auto &tt_run_end) {
    auto lower = std::lower_bound(samples.begin(), samples.end(), tt_run_start);
    samples_idx_sorted.push_back(lower - samples.begin());
    ++n_runs;
  };

  // Cumulative number of BWT-runs per symbols
  auto key_cum = sdsl::cache_file_name(sri::key_trait<t_width>::KEY_BWT_RUN_CUMULATIVE_COUNT, t_config);
  sdsl::int_vector_buffer<> cumulative(key_cum, std::ios::out, buffer_size, log_r);

  auto sigma = psi_rle.sigma();
  for (decltype(sigma) c = 0; c < sigma; ++c) {
    psi_rle.traverse(c, report);
    cumulative.push_back(n_runs);
  }

  samples_idx_sorted.close();
  cumulative.close();
}

template<uint8_t t_width>
void constructSubsamplingBackwardSamplesSortedByAlphabet(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto key_prefix = std::to_string(t_subsample_rate) + "_";
  const std::size_t buffer_size = 1 << 20;

  sdsl::int_vector<> subsamples_idx_to_sorted;
  {
    sdsl::int_vector<> samples_idx_sorted;
    sdsl::load_from_cache(samples_idx_sorted,
                          sri::KeySortedByAlphabet(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX),
                          t_config);

    sdsl::bit_vector subsamples_idx_bv;
    {
      sdsl::int_vector<> subsamples_idx;
      sdsl::load_from_cache(subsamples_idx, key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_config);

      subsamples_idx_bv = sri::constructBitVectorFromIntVector(subsamples_idx, samples_idx_sorted.size(), false);

      subsamples_idx_to_sorted = sdsl::int_vector<>(subsamples_idx.size(), 0, subsamples_idx.width());
    }
    sdsl::bit_vector::rank_1_type rank_subsamples_idx_bv(&subsamples_idx_bv);

    sdsl::int_vector<> subsamples;
    sdsl::load_from_cache(subsamples, key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);

    auto key = sdsl::cache_file_name(
        sri::KeySortedByAlphabet(key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS), t_config);
    sdsl::int_vector_buffer<> subsamples_sorted(key, std::ios::out, buffer_size, subsamples.width());

    key = sdsl::cache_file_name(
        sri::KeySortedByAlphabet(key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX), t_config);
    sdsl::int_vector_buffer<> subsamples_idx_sorted(key, std::ios::out, buffer_size, subsamples_idx_to_sorted.width());

    std::size_t i = 0;
    std::size_t c_subsamples = 0;
    for (auto &&idx: samples_idx_sorted) {
      if (subsamples_idx_bv[idx] == true) {
        auto rnk = rank_subsamples_idx_bv(idx);
        subsamples_sorted.push_back(subsamples[rnk]);
        subsamples_idx_sorted.push_back(i);

        subsamples_idx_to_sorted[rnk] = c_subsamples++;
      }
      ++i;
    }

    subsamples_sorted.close();
    subsamples_idx_sorted.close();
  }

  sdsl::int_vector<> mark_to_sample_idx;
  sdsl::load_from_cache(mark_to_sample_idx,
                        key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
                        t_config);
  auto key = sdsl::cache_file_name(
      sri::KeySortedByAlphabet(key_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX),
      t_config);
  sdsl::int_vector_buffer<> mark_to_sample_idx_sorted(key, std::ios::out, buffer_size, mark_to_sample_idx.width());
  for (const auto &idx: mark_to_sample_idx) {
    mark_to_sample_idx_sorted.push_back(subsamples_idx_to_sorted[idx]);
  }

  mark_to_sample_idx_sorted.close();
}

template<uint8_t t_width>
void constructSubsamplingBackwardSamplesPosition(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  sdsl::int_vector<> samples_pos; // BWT-run starts positions in SA
  sdsl::load_from_cache(samples_pos, sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);

  sdsl::int_vector<> subsamples_idx; // Sub-samples indices
  sdsl::load_from_cache(subsamples_idx, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_config);

  // Compute sub-samples positions
  sdsl::int_vector<> subsamples_pos(subsamples_idx.size(), 0, samples_pos.width());
  std::transform(subsamples_idx.begin(), subsamples_idx.end(), subsamples_pos.begin(),
                 [&samples_pos](auto tt_i) { return samples_pos[tt_i]; });

  sdsl::store_to_cache(subsamples_pos, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);
}

template<typename TGetNextMark, typename TGetNextSubmark, typename TReport>
void computeSubmarksValidity(std::size_t t_r_prime,
                             TGetNextMark t_get_next_mark,
                             TGetNextSubmark t_get_next_submark,
                             TReport t_report) {
  // Marks sampled BWT run heads indices in text and if they are trustworthy
  std::size_t mark = t_get_next_mark();
  std::size_t submark = t_get_next_submark();
  for (std::size_t i = 0, j = 0; i < t_r_prime - 1; ++i) {
    std::size_t next_mark = t_get_next_mark();
    std::size_t next_submark = t_get_next_submark();
    if (next_mark != next_submark) {
      // Report current submark as invalid, and what is the next mark to compute valid area
      t_report(i, submark, next_mark);

      do {
        next_mark = t_get_next_mark();
      } while (next_mark != next_submark);
    }

    submark = next_submark;
  }
}

template<uint8_t t_width>
void constructSubsamplingBackwardMarksValidity(std::size_t t_subsample_rate, sdsl::cache_config &t_config) {
  auto prefix = std::to_string(t_subsample_rate) + "_";

  sdsl::int_vector<> marks;
  sdsl::load_from_cache(marks, sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_config);
  sdsl::int_vector<> sorted_marks_idx;
  sdsl::load_from_cache(sorted_marks_idx, sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX, t_config);
  auto it_marks_idx = sorted_marks_idx.end();
  auto get_next_mark = [&marks, &it_marks_idx]() { return marks[*(--it_marks_idx)]; };

  sdsl::int_vector<> submarks;
  sdsl::load_from_cache(submarks, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config);
  std::sort(submarks.begin(), submarks.end());
  auto it_submarks = submarks.end();
  auto get_next_submark = [&it_submarks]() { return *(--it_submarks); };

  auto r_prime = submarks.size();
  std::vector<std::pair<std::size_t, std::size_t>> validity;
  validity.reserve(r_prime / 4);
  std::size_t max_valid_area = 0;
  auto report = [r_prime, &validity, &max_valid_area](auto tt_i, auto tt_submark, auto tt_next_mark) {
    auto valid_area = tt_submark - tt_next_mark + 1;
    validity.emplace_back(r_prime - tt_i - 1, valid_area);

    if (max_valid_area < valid_area) max_valid_area = valid_area;
  };

  computeSubmarksValidity(r_prime, get_next_mark, get_next_submark, report);

  const std::size_t buffer_size = 1 << 20;
  sdsl::int_vector_buffer<> valid_submarks(
      sdsl::cache_file_name(prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK, t_config),
      std::ios::out,
      buffer_size,
      sdsl::bits::hi(r_prime) + 1);

  sdsl::int_vector_buffer<> valid_areas(
      sdsl::cache_file_name(prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA, t_config),
      std::ios::out,
      buffer_size,
      sdsl::bits::hi(max_valid_area) + 1);

  for (auto it = validity.rbegin(); it != validity.rend(); ++it) {
    auto[idx, area] = *it;
    valid_submarks.push_back(it->first);
    valid_areas.push_back(it->second);
  }

  valid_submarks.close();
  valid_areas.close();
}

#endif //SRI_BENCHMARK_SR_CSA_SR_CSA_H_
