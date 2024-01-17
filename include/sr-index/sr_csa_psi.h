//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 12/12/2023.
//

#ifndef SRI_SR_CSA_PSI_H_
#define SRI_SR_CSA_PSI_H_

#include <cstdint>

#include "r_csa.h"
#include "config.h"
#include "construct.h"
#include "sampling.h"

namespace sri {
template<typename TStorage = GenericStorage,
  typename TAlphabet = Alphabet<>,
  typename TPsiRLE = PsiCoreRLE<>,
  typename TBvMark = sdsl::sd_vector<>,
  typename TMarkToSampleIdx = sdsl::int_vector<>,
  typename TSample = sdsl::int_vector<>,
  typename TBvSampleIdx = sdsl::sd_vector<>,
  typename TRunCumulativeCount = sdsl::int_vector<>>
class SrCSAWithPsiRun : public RCSAWithPsiRun<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample> {
public:
  using Alphabet = TAlphabet;
  using BvMarks = TBvMark;
  using Samples = TSample;
  using BvSamplesIdx = TBvSampleIdx;
  using Base = RCSAWithPsiRun<TStorage, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample>;

  SrCSAWithPsiRun(const TStorage& t_storage, std::size_t t_sr) : Base(t_storage, t_sr) {}

  explicit SrCSAWithPsiRun(std::size_t t_sr) : Base(t_sr) {}

  SrCSAWithPsiRun() = default;

  virtual ~SrCSAWithPsiRun() = default;

  using typename Base::size_type;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, const std::string& name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += Base::serialize(out, v, name);

    written_bytes +=
        this->template serializeItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), out, child, "samples_idx");
    written_bytes +=
        this->template serializeRank<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), out, child, "samples_idx_rank");

    // written_bytes +=
    //   this->template serializeItem<TRunCumulativeCount>(key(SrIndexKey::RUN_CUMULATIVE_COUNT), out, child, "run_cumulative_count");

    return written_bytes;
  }

protected:
  using Base::key;
  void setupKeyNames(const JSON& t_keys) override {
    if (!this->keys_.empty()) return;

    Base::setupKeyNames();
    this->keys_.resize(6);
    using namespace sri::conf;
    // key(SrIndexKey::ALPHABET) = t_keys[kAlphabet];
    // key(SrIndexKey::NAVIGATE) = t_keys[kPsi][kBase];
    key(SrIndexKey::SAMPLES) = key_prefix_ + t_keys[kPsi][kHead][kTextPos].get<std::string>();
    key(SrIndexKey::MARKS) = key_prefix_ + t_keys[kPsi][kTail][kTextPos].get<std::string>();
    key(SrIndexKey::MARK_TO_SAMPLE) = key_prefix_ + t_keys[kPsi][kTail][kTextPosAsc][kLink].get<std::string>();
    key(SrIndexKey::SAMPLES_IDX) = key_prefix_ + t_keys[kPsi][kHead][kIdx].get<std::string>();

    // key(SrIndexKey::RUN_CUMULATIVE_COUNT) = conf::KEY_BWT_RUN_CUMULATIVE_COUNT;
  }

  using typename Base::TSource;

  void loadAllItems(TSource& t_source) override {
    Base::loadAllItems(t_source);

    this->template loadItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    this->template loadBVRank<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);

    // this->template loadItem<TRunCumulativeCount>(key(SrIndexKey::RUN_CUMULATIVE_COUNT), t_source);
  }

  void constructIndex(TSource& t_source) override {
    constructIndex(t_source,
                   [this](auto& tt_source) {
                     return this->constructPhiForRange(tt_source,
                                                       constructGetSampleForRun(tt_source),
                                                       constructSplitRangeInBWTRuns(tt_source),
                                                       constructSplitRunInBWTRuns(tt_source),
                                                       constructUpdateRun(),
                                                       constructIsRunEmpty());
                   });
  }

  using typename Base::Range;
  using typename Base::RangeLF;
  template<typename TPhiRange>
  void constructIndex(TSource& t_source, const TPhiRange& t_phi_range) {
    this->index_.reset(
      new RIndexBase{
        this->constructLF(t_source),
        this->constructComputeDataBackwardSearchStep(
          [](const Range& tt_range, Char tt_c, const RangeLF& tt_next_range, std::size_t tt_step) {
            const auto& run = tt_next_range.start.run;
            return DataBackwardSearchStep{
              tt_step,
              RunDataExt(run.start, tt_c, run.rank, tt_range.start != run.start)
            };
          }
        ),
        this->constructComputeSAValues(
          t_phi_range(t_source),
          this->constructComputeToehold(t_source,
                                        constructGetSampleForRunData(t_source),
                                        constructPsiForRunData(t_source))
        ),
        this->n_,
        [](const auto& tt_step) { return DataBackwardSearchStep{0, RunDataExt()}; },
        this->constructGetSymbol(t_source),
        [](auto tt_seq_size) { return Range{0, tt_seq_size}; },
        this->constructIsRangeEmpty()
      }
    );
  }

  using typename Base::RunData;
  using typename Base::Char;
  struct RunDataExt : RunData {
    Char c;
    std::size_t partial_rank;
    bool is_run_start;

    explicit RunDataExt(std::size_t t_pos = 0,
                        Char t_c = 0,
                        std::size_t t_partial_rank = 0,
                        bool t_is_run_start = false)
      : RunData(t_pos), c{t_c}, partial_rank{t_partial_rank}, is_run_start{t_is_run_start} {}
  };

  using DataBackwardSearchStep = sri::DataBackwardSearchStep<RunDataExt>;

  using typename Base::Value;
  auto constructGetSample(TSource& t_source) {
    auto cref_run_cum_c = this->template loadItem<TRunCumulativeCount>(key(SrIndexKey::RUN_CUMULATIVE_COUNT), t_source);
    auto cref_bv_sample_idx = this->template loadItem<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto bv_sample_idx_rank = this->template loadBVRank<TBvSampleIdx>(key(SrIndexKey::SAMPLES_IDX), t_source, true);
    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);

    auto get_sample = [cref_run_cum_c, cref_bv_sample_idx, bv_sample_idx_rank, cref_samples](
      auto tt_char,
      auto tt_partial_rank,
      bool tt_is_run_start
    ) -> std::optional<Value> {
      if (!tt_is_run_start) return std::nullopt;

      std::size_t cc = 0 < tt_char ? cref_run_cum_c.get()[tt_char - 1] : 0;
      auto idx = cc + tt_partial_rank;
      if (!cref_bv_sample_idx.get()[idx]) return std::nullopt;

      auto rank = bv_sample_idx_rank(idx);
      return cref_samples.get()[rank];
    };

    return get_sample;
  }

  auto constructGetSampleForRunData(TSource& t_source) {
    auto get_sample = constructGetSample(t_source);
    return [get_sample](const RunDataExt& tt_run_data) {
      return get_sample(tt_run_data.c, tt_run_data.partial_rank, tt_run_data.is_run_start);
    };
  }

  struct Run {
    std::size_t start = 0;
    std::size_t end = 0;
    Char c = 0;
    std::size_t partial_rank = 0;
    bool is_run_start = false;

    bool operator<(const Run& rhs) const {
      return start < rhs.start;
    }
  };

  auto constructUpdateRun() {
    return [](Run& tt_run) {
      ++tt_run.start;
      tt_run.is_run_start = false;
    };
  }

  auto constructIsRunEmpty() {
    return [](const Run& tt_run) { return !(tt_run.start < tt_run.end); };
  }

  auto constructGetSampleForRun(TSource& t_source) {
    auto get_sample = constructGetSample(t_source);

    return [get_sample](const Run& tt_run) {
      return get_sample(tt_run.c, tt_run.partial_rank, tt_run.is_run_start);
    };
  }

  auto constructCreateRun() {
    auto create_run = [](auto tt_first, auto tt_last, auto tt_c, auto tt_partial_rank, auto tt_is_first) {
      return Run{tt_first, tt_last, tt_c, tt_partial_rank, tt_is_first};
    };

    return create_run;
  }

  auto constructSplitRangeInBWTRuns(TSource& t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(key(SrIndexKey::NAVIGATE), t_source, true);
    auto create_run = constructCreateRun();

    return [cref_psi_core, create_run](const auto& tt_range) -> auto {
      const auto& [first, last] = tt_range;
      return cref_psi_core.get().splitInSortedRuns(first, last, create_run);
    };
  }

  auto constructSplitRunInBWTRuns(TSource& t_source) {
    return Base::constructSplitRunInBWTRuns(t_source, constructCreateRun());
  }

  auto constructPsiForRunData(TSource& t_source) {
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
    auto get_c = [cref_alphabet](auto tt_index) { return computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));

    auto psi = Psi(psi_select, get_c, cumulative);

    return [psi](const RunDataExt& tt_run_data) { return psi(tt_run_data.pos); };
  }

  template<typename TGetSampleRun, typename TSplitRangeInBWTRuns, typename TSplitRunInBWTRuns, typename TUpdateRun,
    typename TIsRunEmpty>
  auto constructPhiForRange(TSource& t_source,
                            const TGetSampleRun& t_get_sample,
                            const TSplitRangeInBWTRuns& t_split_range,
                            const TSplitRunInBWTRuns& t_split_run,
                            const TUpdateRun& t_update_run,
                            const TIsRunEmpty& t_is_run_empty) {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(key(SrIndexKey::MARKS), t_source, true);
    auto successor = CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(key(SrIndexKey::MARK_TO_SAMPLE), t_source);
    auto get_mark_to_sample_idx = RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, false);

    auto cref_samples = this->template loadItem<TSample>(key(SrIndexKey::SAMPLES), t_source);
    auto get_sample = RandomAccessForCRefContainer(cref_samples);
    SampleValidatorDefault sample_validator_default;

    auto phi = buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);
    auto phi_simple = [phi](const auto& tt_prev_value) { return phi(tt_prev_value).first; };

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

  std::size_t subsample_rate_ = 1;
  std::string key_prefix_;
};

template<uint8_t t_width, typename TSamples, typename TBvMark>
void constructBaseSrCSAWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

template<typename... TArgs>
void construct(SrCSAWithPsiRun<TArgs...>& t_index, const std::string& t_data_path, Config& t_config) {
  using Index = SrCSAWithPsiRun<TArgs...>;
  using width = typename Index::Alphabet::int_width;
  using namespace conf;
  const auto& keys = t_config.keys;

  constructRCSAWithPsiRuns<width, typename Index::BvMarks>(t_data_path, t_config);

  auto subsample_rate = t_index.SubsampleRate();

  constructBaseSrCSAWithPsiRuns<width, typename Index::Samples, typename Index::BvMarks>(subsample_rate, t_config);

  auto n = sizeIntVector<width>(t_config, keys[kBWT][kBase]);

  auto prefix = std::to_string(subsample_rate) + "_";

  if (
    auto key = prefix + keys[kPsi][kHead][kIdx].get<std::string>();
    !sdsl::cache_file_exists<typename Index::BvSamplesIdx>(key, t_config)
  ) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructBitVectorFromIntVector<typename Index::BvSamplesIdx>(key, t_config, n, false);
  }

  t_index.load(t_config);
}

template<typename TSamples>
void constructSubsamplesForPhiForwardWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

auto constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

void constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

template<uint8_t t_width, typename TSamples = sdsl::int_vector<>, typename TBvMark = sdsl::sd_vector<>>
void constructBaseSrCSAWithPsiRuns(const std::size_t t_subsample_rate, Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;

  auto n = sizeIntVector<t_width>(t_config, keys[kBWT][kBase]);

  const auto prefix = std::to_string(t_subsample_rate) + "_";

  // Sort samples (Psi-run head) by its text positions
  if (!cache_file_exists(keys[kPsi][kHead][kTextPosAsc][kIdx], t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSortedIndices(keys[kPsi][kHead][kTextPos], t_config, keys[kPsi][kHead][kTextPosAsc][kIdx]);
  }

  // Construct subsampling backward of samples (text positions of Psi-run first letter)
  if (!sdsl::cache_file_exists<TSamples>(prefix + keys[kPsi][kHead][kTextPos].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplesForPhiForwardWithPsiRuns<TSamples>(t_subsample_rate, t_config);
  }

  // Construct subsampling backward of marks (text positions of Psi-run last letter)
  if (!cache_file_exists(prefix + keys[kPsi][kTail][kTextPosAsc][kLink].get<std::string>(), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(t_subsample_rate, t_config);
  }

  // Construct successor on the text positions of sub-sampled Psi-run last letter
  if (
    auto key = prefix + keys[kPsi][kTail][kTextPos].get<std::string>();
    !sdsl::cache_file_exists<TBvMark>(key, t_config)
  ) {
    auto event = sdsl::memory_monitor::event("Successor");
    constructBitVectorFromIntVector<TBvMark>(key, t_config, n, false);
  }
}

inline auto constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(const std::size_t t_subsample_rate,
                                                                        Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;
  const auto prefix = std::to_string(t_subsample_rate) + "_";

  // Samples
  sdsl::int_vector<> samples; // Psi-run starts positions in text
  sdsl::load_from_cache(samples, keys[kPsi][kHead][kTextPos], t_config);

  sdsl::int_vector<> sorted_samples_idx;
  sdsl::load_from_cache(sorted_samples_idx, keys[kPsi][kHead][kTextPosAsc][kIdx], t_config);

  // We must sub-sample the samples associated to the first and last marks in the text
  const auto req_samples_idx = getExtremes(t_config, keys[kPsi][kTail][kTextPosAsc][kLink]);

  // Compute sub-sampled indices for sampled values (sorted by Psi positions)
  const auto subsamples_idx = computeSortedSampling(t_subsample_rate,
                                                    std::reverse_iterator(sorted_samples_idx.end()),
                                                    // sorted_samples_idx.begin(),
                                                    std::reverse_iterator(sorted_samples_idx.begin()),
                                                    // sorted_samples_idx.end(),
                                                    samples,
                                                    req_samples_idx);
  sdsl::store_to_cache(subsamples_idx, prefix + keys[kPsi][kHead][kIdx].get<std::string>(), t_config);

  // Compute sub-samples
  sdsl::int_vector<> subsamples(subsamples_idx.size(), 0, samples.width());
  std::transform(subsamples_idx.begin(),
                 subsamples_idx.end(),
                 subsamples.begin(),
                 [&samples](auto tt_i) { return samples[tt_i]; });
  sri::store_to_cache(subsamples, prefix + keys[kPsi][kHead][kTextPos].get<std::string>(), t_config, true);

  return subsamples;
}

template<typename TSamples>
void constructSubsamplesForPhiForwardWithPsiRuns(const std::size_t t_subsample_rate, Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;
  const auto prefix = std::to_string(t_subsample_rate) + "_";
  const auto key = prefix + keys[kPsi][kHead][kTextPos].get<std::string>();

  sdsl::int_vector<> subsamples_iv;
  if (!sdsl::cache_file_exists<sdsl::int_vector<>>(key, t_config)) {
    subsamples_iv = constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(t_subsample_rate, t_config);
  } else {
    sdsl::load_from_cache(subsamples_iv, key, t_config, true);
  }

  if (!std::is_same_v<TSamples, sdsl::int_vector<>>) {
    auto subsamples = construct<TSamples>(subsamples_iv);
    sdsl::store_to_cache(subsamples, key, t_config, true);
  }
}

inline auto computeSampleToMarkLinksForPhiForwardWithPsiRuns(const std::string& t_prefix, const Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;

  // Sub-sampled indices of samples
  sdsl::int_vector<> subsamples_idx;
  sdsl::load_from_cache(subsamples_idx, t_prefix + keys[kPsi][kHead][kIdx].get<std::string>(), t_config);

  sdsl::int_vector<> subsample_to_mark_links(subsamples_idx.size(), 0, subsamples_idx.width());

  const auto r = sizeIntVector<0>(t_config, keys[kPsi][kHead][kTextPos]);

  // Compute links from samples to marks
  for (int i = 0; i < subsamples_idx.size(); ++i) {
    subsample_to_mark_links[i] = (subsamples_idx[i] + r - 1) % r;
  }

  return subsample_to_mark_links;
}

inline auto computeSubmarksForPhiForwardWithPsiRuns(const std::string& t_prefix, const Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;

  // Links from sub-samples to their corresponding mark (actually, to rank of mark in BWT)
  const auto subsample_to_mark_links = computeSampleToMarkLinksForPhiForwardWithPsiRuns(t_prefix, t_config);
  const auto r_prime = subsample_to_mark_links.size();

  sdsl::int_vector<> marks;
  sdsl::load_from_cache(marks, keys[kPsi][kTail][kTextPos], t_config);

  auto submarks = sdsl::int_vector(r_prime, 0, marks.width());
  std::transform(subsample_to_mark_links.begin(),
                 subsample_to_mark_links.end(),
                 submarks.begin(),
                 [&marks](auto tt_i) { return marks[tt_i]; });

  return submarks;
}

inline void constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(const std::size_t t_subsample_rate,
                                                                      Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;

  const auto prefix = std::to_string(t_subsample_rate) + "_";

  // Text positions of marks indices associated to sub-samples, i.e., text positions of sub-sampled marks.
  // Note that the submarks are sorted by its associated submark position, not by its positions in Psi
  sdsl::int_vector<> submarks = computeSubmarksForPhiForwardWithPsiRuns(prefix, t_config);

  const auto r_prime = submarks.size();
  const auto log_r_prime = sdsl::bits::hi(r_prime) + 1;

  // Links from sub-sampled marks (sorted by text position) to sub-samples indices.
  // Note that, initially, these are the indices of sub-sample in Psi.
  sdsl::int_vector<> submark_to_subsample_links(r_prime, 0, log_r_prime);
  std::iota(submark_to_subsample_links.begin(), submark_to_subsample_links.end(), 0);

  // Sort indexes by text positions of its marks, becoming in the links from the sub-sampled marks to sub-sampled samples.
  std::sort(submark_to_subsample_links.begin(),
            submark_to_subsample_links.end(),
            [&submarks](const auto& tt_a, const auto& tt_b) { return submarks[tt_a] < submarks[tt_b]; });

  sdsl::store_to_cache(submarks, prefix + keys[kPsi][kTail][kTextPos].get<std::string>(), t_config);
  sdsl::store_to_cache(submark_to_subsample_links,
                       prefix + keys[kPsi][kTail][kTextPosAsc][kLink].get<std::string>(),
                       t_config);
}
}

#endif //SRI_SR_CSA_PSI_H_
