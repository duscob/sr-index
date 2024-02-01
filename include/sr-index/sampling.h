//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/21.
//

#ifndef SRI_SAMPLING_H_
#define SRI_SAMPLING_H_

#include <cstddef>
#include <iterator>
#include <vector>
#include <type_traits>
#include <algorithm>

namespace sri {
template<typename TInputIt, typename TValues, typename TContainer>
auto computeSampling(std::size_t t_s,
                     TInputIt t_first_idx,
                     TInputIt t_last_idx,
                     const TValues& t_values,
                     const TContainer& t_req_idxs) {
  auto n = std::distance(t_first_idx, t_last_idx);
  std::vector<typename TInputIt::value_type> samples;
  samples.reserve(n / 2);

  if (n == 0) return samples;

  samples.emplace_back(*t_first_idx);
  auto last_sampled_value = t_values[*t_first_idx];

  if (n == 1) return samples;

  auto prev_idx = t_first_idx + 1;
  auto prev = t_values[*prev_idx];
  for (auto curr_idx = prev_idx + 1; curr_idx != t_last_idx; ++curr_idx) {
    auto curr = t_values[*curr_idx];

    if (t_s < (curr < last_sampled_value ? last_sampled_value - curr : curr - last_sampled_value)
      || std::find(t_req_idxs.begin(), t_req_idxs.end(), *prev_idx) != t_req_idxs.end()) {
      samples.emplace_back(*prev_idx);
      last_sampled_value = prev;
    }

    prev_idx = curr_idx;
    prev = curr;
  }

  samples.emplace_back(*(prev_idx));

  return samples;
}

template<typename TInputIt, typename TValues, typename TContainer>
auto computeSortedSampling(std::size_t t_s,
                           TInputIt t_first_idx,
                           TInputIt t_last_idx,
                           const TValues& t_values,
                           const TContainer& t_req_idxs) {
  auto subsamples_idx_vec = computeSampling(t_s, t_first_idx, t_last_idx, t_values, t_req_idxs);

  std::sort(subsamples_idx_vec.begin(), subsamples_idx_vec.end());

  auto r = t_values.size();
  auto log_r = sdsl::bits::hi(r) + 1;
  auto subsamples_idx = sdsl::int_vector<>(subsamples_idx_vec.size(), 0, log_r);
  std::copy(subsamples_idx_vec.begin(), subsamples_idx_vec.end(), subsamples_idx.begin());

  return subsamples_idx;
}
}

#endif //SRI_SAMPLING_H_
