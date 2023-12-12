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

template<uint8_t t_width>
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
  if (!cache_file_exists(prefix + to_string(keys[kPsi][kHead][kIdx]), t_config)) {
    auto event = sdsl::memory_monitor::event("Subsampling");
    constructSubsamplingBackwardSamplesForPhiForwardWithPsiRuns(t_subsample_rate, t_config);
  }

  // // Construct subsampling backward of marks (text positions of Psi-run last letter)
  // if (!cache_file_exists(prefix + to_string(keys[kPsi][kTail][kTextPosAsc][kLink]), t_config)) {
  //   auto event = sdsl::memory_monitor::event("Subsampling");
  //   constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns<t_width>(t_subsample_rate, t_config);
  // }
  //
  // // Construct successor on the text positions of sub-sampled BWT-run last letter
  // if (!sdsl::cache_file_exists<TBvMark>(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config)) {
  //   auto event = sdsl::memory_monitor::event("Successor");
  //   constructBitVectorFromIntVector<TBvMark>(prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config, n, false);
  // }
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

// template<uint8_t t_width>
// void constructSubsamplingBackwardMarksForPhiForwardWithPsiRuns(const std::size_t t_subsample_rate, Config& t_config) {
//   const auto prefix = std::to_string(t_subsample_rate) + "_";
//
//   // Text positions of marks indices associated to sub-samples, i.e., text positions of sub-sampled marks
//   sdsl::int_vector<> subsampled_mark_text_pos;
//   std::size_t r_prime;
//
//   {
//     // Compute sub-sample to mark links
//
//     // Links from sub-samples to their corresponding mark (actually, to rank of mark in BWT)
//     auto subsample_to_mark_links = computeSampleToMarkLinksForPhiForward<t_width>(prefix, t_config);
//     r_prime = subsample_to_mark_links.size();
//
//     sdsl::int_vector<> bwt_run_ends_text_pos;
//     sdsl::load_from_cache(bwt_run_ends_text_pos, conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config);
//
//     subsampled_mark_text_pos = sdsl::int_vector(r_prime, 0, bwt_run_ends_text_pos.width());
//     std::transform(subsample_to_mark_links.begin(),
//                    subsample_to_mark_links.end(),
//                    subsampled_mark_text_pos.begin(),
//                    [&bwt_run_ends_text_pos](auto tt_i) { return bwt_run_ends_text_pos[tt_i]; });
//   }
//
//   auto log_r_prime = sdsl::bits::hi(r_prime) + 1;
//
//   // Links from sub-sampled marks (sorted by text position) to sub-samples indices. Note that, initially, these are the indices of sub-sample in BWT.
//   sdsl::int_vector<> subsampled_mark_to_subsample_links(r_prime, 0, log_r_prime);
//   std::iota(subsampled_mark_to_subsample_links.begin(), subsampled_mark_to_subsample_links.end(), 0);
//
//   // Sort indexes by text positions of its marks, becoming in the links from the sub-sampled marks to sub-sampled samples.
//   std::sort(subsampled_mark_to_subsample_links.begin(),
//             subsampled_mark_to_subsample_links.end(),
//             [&subsampled_mark_text_pos](const auto &tt_a, const auto &tt_b) {
//               return subsampled_mark_text_pos[tt_a] < subsampled_mark_text_pos[tt_b];
//             });
//
//   sdsl::store_to_cache(subsampled_mark_text_pos, prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST, t_config);
//
//   sdsl::store_to_cache(subsampled_mark_to_subsample_links,
//                        prefix + conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
//                        t_config);
// }
}

#endif //SRI_SR_CSA_PSI_H_
