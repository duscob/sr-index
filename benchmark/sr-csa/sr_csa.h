//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_SR_CSA_H_
#define SRI_BENCHMARK_SR_CSA_SR_CSA_H_

#include <optional>

#include "sr-index/sampling.h"

#include "csa.h"

struct RunSimpleSrCSA {
  std::size_t first;
  std::size_t last;
};

bool operator<(const RunSimpleSrCSA &lhs, const RunSimpleSrCSA &rhs) {
  return lhs.first < rhs.first;
}

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSampleIdx = sdsl::sd_vector<>,
    typename TRun = RunSimpleSrCSA>
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
    written_bytes += this->template serializeItem<TAlphabet>(
        sri::key_trait<t_width>::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += this->template serializeItem<TBvMark>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks_select");

    written_bytes += this->template serializeItem<TMarkToSampleIdx>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
        out,
        child,
        "mark_to_sample");

    written_bytes += this->template serializeItem<TSample>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, out, child, "samples");

    written_bytes += this->template serializeItem<TBvSampleIdx>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, out, child, "samples_idx");

    return written_bytes;
  }

 protected:

  using typename BaseClass::TSource;

  using typename BaseClass::Value;
  using TFnGetSample = std::function<std::optional<Value>(std::size_t)>;
  virtual TFnGetSample constructGetSampleForSAIdx(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);
    auto get_run = [cref_psi_core](auto tt_sa_pos) {
      auto[run, run_start] = cref_psi_core.get().rankSoftRun(tt_sa_pos);
      return std::make_pair(run - 1, run_start);
    };

    auto cref_bv_sample_idx = this->template loadItem<TBvSampleIdx>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_source, true);
    auto is_run_sampled = [cref_bv_sample_idx](auto tt_i) {
      return cref_bv_sample_idx.get()[tt_i];
    };

    auto cref_samples = this->template loadItem<TSample>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto bv_sample_idx_rank = this->template loadBVRank<TBvSampleIdx>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_IDX, t_source, true);
    auto get_sample = [cref_samples, bv_sample_idx_rank](auto tt_run) {
      auto sampled_run = bv_sample_idx_rank(tt_run);
      return cref_samples.get()[sampled_run];
    };

    return sri::GetSampleForSAPosition(get_run, is_run_sampled, get_sample);
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
    auto cref_psi_core = this->template loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);
    auto create_run = constructCreateRun();

    auto split = [cref_psi_core, create_run](const auto &tt_range) -> auto {
      const auto &[first, last] = tt_range;
      return cref_psi_core.get().splitInSortedRuns(first, last, create_run);
    };

    return split;
  }

  using TFnSplitRun = std::function<Runs(TRun)>;
  virtual TFnSplitRun constructSplitRunInBWTRuns(TSource &t_source) {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);
    auto cref_alphabet = this->template loadItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, t_source);
    auto create_run = constructCreateRun();

    auto split = [cref_psi_core, cref_alphabet, create_run](const auto &tt_run) -> auto {
      const auto &[first, last] = tt_run;
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

  using typename BaseClass::TFnPhiForRange;
  TFnPhiForRange constructPhiForRange(TSource &t_source) override {
    auto bv_mark_rank = this->template loadBVRank<TBvMark>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_source, true);
    auto bv_mark_select = this->template loadBVSelect<TBvMark>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_source, true);
    auto successor = sri::CircularSoftSuccessor(bv_mark_rank, bv_mark_select, this->n_);

    auto cref_mark_to_sample_idx = this->template loadItem<TMarkToSampleIdx>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_source);
    auto get_mark_to_sample_idx = sri::RandomAccessForTwoContainersDefault(cref_mark_to_sample_idx, false);

    auto cref_samples = this->template loadItem<TSample>(
        key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto get_sample = sri::RandomAccessForCRefContainer(cref_samples);
    sri::SampleValidatorDefault sample_validator_default;

    auto phi = sri::buildPhiForward(successor, get_mark_to_sample_idx, get_sample, sample_validator_default, this->n_);

    auto split_range = constructSplitRangeInBWTRuns(t_source);
    auto split_run = constructSplitRunInBWTRuns(t_source);
    auto get_sample_x_sa_idx = constructGetSampleForSAIdx(t_source);

    return sri::buildPhiForwardForRangeSimple(
        phi, split_range, split_run, get_sample_x_sa_idx, subsample_rate_, this->n_);
  }

  using typename BaseClass::TFnComputeToehold;
  TFnComputeToehold constructComputeToeholdForPhiForward(TSource &t_source) override {
    auto cref_psi_core = this->template loadItem<TPsiRLE>(sdsl::conf::KEY_PSI, t_source, true);

    auto get_sample = constructGetSampleForSAIdx(t_source);

    auto psi_select = [cref_psi_core](auto tt_c, auto tt_rnk) { return cref_psi_core.get().select(tt_c, tt_rnk); };
    auto cref_alphabet = this->template loadItem<TAlphabet>(sri::key_trait<t_width>::KEY_ALPHABET, t_source);
    auto get_c = [cref_alphabet](auto tt_index) { return sri::computeCForSAIndex(cref_alphabet.get().C, tt_index); };
    auto cumulative = sri::RandomAccessForCRefContainer(std::cref(cref_alphabet.get().C));
    auto psi = sri::Psi(psi_select, get_c, cumulative);

    auto get_sa_value_for_bwt_run_start = sri::buildComputeSAValueForward(get_sample, psi, this->n_);

    return sri::buildComputeToeholdForPhiForward(cref_psi_core, get_sa_value_for_bwt_run_start);
  }

  std::size_t subsample_rate_ = 1;
  std::string key_prefix_;
};

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBvMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>,
    typename TBvSamplePos = sdsl::sd_vector<>,
    typename TRun = RunSimpleSrCSA>
class SrCSAWithBv : public SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos, TRun> {
 public:
  using BaseClass = SrCSA<t_width, TAlphabet, TPsiRLE, TBvMark, TMarkToSampleIdx, TSample, TBvSamplePos, TRun>;

  SrCSAWithBv(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr) : BaseClass(t_storage, t_sr) {}

  using typename BaseClass::size_type;

  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, const std::string &name) const override {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += this->template serializeItem<TAlphabet>(
        sri::key_trait<t_width>::KEY_ALPHABET, out, child, "alphabet");

    written_bytes += this->template serializeItem<TPsiRLE>(sdsl::conf::KEY_PSI, out, child, "psi");

    written_bytes += this->template serializeItem<TBvMark>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks");
    written_bytes += this->template serializeRank<TBvMark>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks_rank");
    written_bytes += this->template serializeSelect<TBvMark>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, out, child, "marks_select");

    written_bytes += this->template serializeItem<TMarkToSampleIdx>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
        out,
        child,
        "mark_to_sample");

    written_bytes += this->template serializeItem<TSample>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, out, child, "samples");

    written_bytes += this->template serializeItem<TBvSamplePos>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, out, child, "samples_pos");

    return written_bytes;
  }

 protected:

  using typename BaseClass::TSource;

  using typename BaseClass::TFnGetSample;
  TFnGetSample constructGetSampleForSAIdx(TSource &t_source) override {
    auto cref_samples = this->template loadItem<TSample>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_source);
    auto cref_bv_sample_pos = this->template loadItem<TBvSamplePos>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_source, true);
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_source, true);

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
    auto bv_sample_pos_rank = this->template loadBVRank<TBvSamplePos>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_source, true);
    auto bv_sample_pos_select = this->template loadBVSelect<TBvSamplePos>(
        this->key_prefix_ + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST, t_source, true);

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

template<uint8_t t_width, typename TBVMark>
void constructSrCSACommons(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBVMark, typename TMarkToSampleIdx, typename TSample, typename TBVOriginalSampleIdx>
void construct(SrCSA<t_width, TAlphabet, TPsiCore, TBVMark, TMarkToSampleIdx, TSample, TBVOriginalSampleIdx> &t_index,
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
    if (!sdsl::cache_file_exists<TBVOriginalSampleIdx>(key, t_config)) {
      sri::constructBitVectorFromIntVector<TBVOriginalSampleIdx>(key, t_config, r);
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
      sri::constructBitVectorFromIntVector<TBvSamplePos>(key, t_config, n);
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
      sri::constructBitVectorFromIntVector<TBvMark>(key, t_config, n);
    }
  }
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

#endif //SRI_BENCHMARK_SR_CSA_SR_CSA_H_
