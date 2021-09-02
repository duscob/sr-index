//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/30/21.
//

#ifndef SRI_BENCHMARK_SR_CSA_SR_CSA_H_
#define SRI_BENCHMARK_SR_CSA_SR_CSA_H_

#include "sr-index/sampling.h"

#include "csa.h"

template<uint8_t t_width = 8,
    typename TAlphabet = sdsl::byte_alphabet,
    typename TPsiRLE = sri::PsiCoreRLE<>,
    typename TBVMark = sdsl::sd_vector<>,
    typename TMarkToSampleIdx = sdsl::int_vector<>,
    typename TSample = sdsl::int_vector<>>
class SrCSA : public CSA<t_width, TAlphabet, TPsiRLE, TBVMark, TMarkToSampleIdx, TSample> {

 public:

  explicit SrCSA(std::reference_wrapper<ExternalStorage> t_storage, std::size_t t_sr)
      : CSA<t_width, TAlphabet, TPsiRLE, TBVMark, TMarkToSampleIdx, TSample>(t_storage), subsample_rate_{t_sr} {
  }

  void load(sdsl::cache_config t_config) override {
  }

  auto SubsampleRate() const { return subsample_rate_; }

 protected:

  std::size_t subsample_rate_ = 1;
};

template<uint8_t t_width>
void constructSubsamplingBackwardSamples(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width>
void constructSubsamplingBackwardMarks(std::size_t t_subsample_rate, sdsl::cache_config &t_config);

template<uint8_t t_width, typename TAlphabet, typename TPsiCore, typename TBVMark, typename TMarkToSampleIdx, typename TSample>
void construct(SrCSA<t_width, TAlphabet, TPsiCore, TBVMark, TMarkToSampleIdx, TSample> &t_index,
               sdsl::cache_config &t_config) {
  std::size_t n;
  {
    sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
    n = bwt_buf.size();
  }

  auto subsample_rate = t_index.SubsampleRate();
  auto prefix = std::to_string(subsample_rate) + "_";

  {
    // Sort samples (BWT-run last letter) by its text positions
    auto event = sdsl::memory_monitor::event("Subsampling");
    const auto key = sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX;
    if (!sdsl::cache_file_exists<TBVMark>(key, t_config)) {
      sri::constructSortedIndices(sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config, key);
    }
  }

  {
    // Construct subsampling backward of samples (text positions of BWT-run last letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    if (!sdsl::cache_file_exists(prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_SAMPLED, t_config)) {
      constructSubsamplingBackwardSamples<t_width>(subsample_rate, t_config);
    }
  }

  {
    // Construct subsampling backward of marks (text positions of BWT-run first letter)
    auto event = sdsl::memory_monitor::event("Subsampling");
    if (!sdsl::cache_file_exists(
        prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX_SAMPLED, t_config)) {
      constructSubsamplingBackwardMarks<t_width>(subsample_rate, t_config);
    }
  }

  {
    // Construct successor on the text positions of sub-sampled BWT-run last letter
    auto event = sdsl::memory_monitor::event("Successor");
    const auto key = prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST_SAMPLED;
    if (!sdsl::cache_file_exists<TBVMark>(key, t_config)) {
      sri::constructBitVectorFromIntVector<TBVMark>(key, t_config, n);
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
    sdsl::store_to_cache(subsamples_idx, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_SAMPLED, t_config);
  }

  {
    // Compute sub-samples
    sdsl::int_vector<> subsamples(subsamples_idx.size(), 0, samples.width());
    std::transform(subsamples_idx.begin(), subsamples_idx.end(), subsamples.begin(),
                   [&samples](auto tt_i) { return samples[tt_i]; });

    sdsl::store_to_cache(subsamples, prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS_SAMPLED, t_config);
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
                       prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST_SAMPLED,
                       t_config);

  sdsl::store_to_cache(subsampled_mark_to_subsample_links,
                       prefix + sri::key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX_SAMPLED,
                       t_config);
}

template<uint8_t t_width>
auto computeSampleToMarkLinks(const std::string &t_prefix, sdsl::cache_config &t_config) {
  // Sub-sampled indices of samples
  sdsl::int_vector<> subsamples_idx;
  sdsl::load_from_cache(subsamples_idx, t_prefix + sri::key_trait<t_width>::KEY_BWT_RUN_FIRST_SAMPLED, t_config);

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

#endif //SRI_BENCHMARK_SR_CSA_SR_CSA_H_
