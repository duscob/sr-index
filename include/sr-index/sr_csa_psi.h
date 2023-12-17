//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 12/12/2023.
//

#ifndef SRI_SR_CSA_PSI_H_
#define SRI_SR_CSA_PSI_H_

#include <cstdint>

#include "config.h"
#include "construct.h"
#include "sampling.h"

namespace sri {
void constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

void constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(std::size_t t_subsample_rate, Config& t_config);

template<uint8_t t_width, typename TBvMark>
void constructSrCSACommonsWithPsiRuns(const std::size_t t_subsample_rate, Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;
  auto n = sizeIntVector<t_width>(t_config, t_config.keys[kBWT][kBase]);

  const auto prefix = std::to_string(t_subsample_rate) + "_";

  // Sort samples (Psi-run head) by its text positions
  if (!cache_file_exists(keys[kPsi][kHead][kTextPosAsc][kIdx], t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSortedIndices(keys[kPsi][kHead][kTextPos], t_config, keys[kPsi][kHead][kTextPosAsc][kIdx]);
  }

  // Construct subsampling backward of samples (text positions of Psi-run first letter)
  if (!cache_file_exists(prefix + to_string(keys[kPsi][kHead][kTextPos]), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(t_subsample_rate, t_config);
  }

  // Construct subsampling backward of marks (text positions of Psi-run last letter)
  if (!cache_file_exists(prefix + to_string(keys[kPsi][kTail][kTextPosAsc][kLink]), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(t_subsample_rate, t_config);
  }

  // Construct successor on the text positions of sub-sampled Psi-run last letter
  if (!sdsl::cache_file_exists<TBvMark>(prefix + to_string(keys[kPsi][kTail][kTextPos]), t_config)) {
    auto event = sdsl::memory_monitor::event("Successor");
    constructBitVectorFromIntVector<TBvMark>(prefix + to_string(keys[kPsi][kTail][kTextPos]), t_config, n, false);
  }
}

inline void constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(const std::size_t t_subsample_rate,
                                                                        Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;
  const auto prefix = std::to_string(t_subsample_rate) + "_";

  // Samples
  sdsl::int_vector<> samples; // BWT-run starts positions in text
  sdsl::load_from_cache(samples, keys[kPsi][kHead][kTextPos], t_config);

  sdsl::int_vector<> sorted_samples_idx;
  sdsl::load_from_cache(sorted_samples_idx, keys[kPsi][kHead][kTextPosAsc][kIdx], t_config);

  // TODO Check if it is required to sampled the extremes
  // We must sub-sample the samples associated to the first and last marks in the text
  // const auto req_samples_idx = getExtremes(t_config, keys[kPsi][kTail][kTextPosAsc][kLink]);
  constexpr std::array<std::size_t, 2> req_samples_idx{};

  // Compute sub-sampled indices for sampled values (sorted by BWT positions)
  const auto subsamples_idx = computeSortedSampling(t_subsample_rate,
                                                    // std::reverse_iterator(sorted_samples_idx.end()),
                                                    sorted_samples_idx.begin(),
                                                    // std::reverse_iterator(sorted_samples_idx.begin()),
                                                    sorted_samples_idx.end(),
                                                    samples,
                                                    req_samples_idx);
  sdsl::store_to_cache(subsamples_idx, prefix + to_string(keys[kPsi][kHead][kIdx]), t_config);

  // Compute sub-samples
  sdsl::int_vector<> subsamples(subsamples_idx.size(), 0, samples.width());
  std::transform(subsamples_idx.begin(),
                 subsamples_idx.end(),
                 subsamples.begin(),
                 [&samples](auto tt_i) { return samples[tt_i]; });
  sdsl::store_to_cache(subsamples, prefix + to_string(keys[kPsi][kHead][kTextPos]), t_config);
}

inline auto computeSampleToMarkLinksForPhiForwardWithPsiRuns(const std::string& t_prefix, const Config& t_config) {
  using namespace conf;
  const auto& keys = t_config.keys;

  // Sub-sampled indices of samples
  sdsl::int_vector<> subsamples_idx;
  sdsl::load_from_cache(subsamples_idx, t_prefix + to_string(keys[kPsi][kHead][kIdx]), t_config);

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

  // Text positions of marks indices associated to sub-samples, i.e., text positions of sub-sampled marks
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

  sdsl::store_to_cache(submarks, prefix + to_string(keys[kPsi][kTail][kTextPos]), t_config);
  sdsl::store_to_cache(submark_to_subsample_links, prefix + to_string(keys[kPsi][kTail][kTextPosAsc][kLink]), t_config);
}
}

#endif //SRI_SR_CSA_PSI_H_
