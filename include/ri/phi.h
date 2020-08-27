//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_PHI_H_
#define RI_PHI_H_

#include <cstddef>
#include <cassert>
#include <iterator>
#include <experimental/optional>

namespace ri {

template<typename TPartitioner, typename TPred, typename TPredToRun, typename TSampledTail>
class Phi {
 public:
/*
 * Phi function. Phi(SA[0]) is undefined
 */
  std::size_t operator()(std::size_t prev_value) {

    assert(prev_value != bwt_size_ - 1);

    auto p = pred_(prev_value);

    //pred_pos is the rank of the predecessor of prev_value (circular)
    auto pred_pos = p.first;
//    std::size_t pred_pos = pred.predecessor_rank_circular(prev_value);
//    assert(pred_pos <= r - 1);

    //the actual predecessor
    auto pred_value = p.second;
//    ulint pred_value = pred.select(pred_pos);
//    assert(pred_pos < r - 1 or pred_value == bwt_size_ - 1);

    //distance from predecessor
    auto delta = pred_value < prev_value ? prev_value - pred_value : prev_value + 1;

    //cannot fall on first run: this can happen only if I call Phi(SA[0])
    assert(0 < pred_to_run_(pred_pos));
    auto run = pred_to_run_(pred_pos) - 1;

    //sample at the end of previous run
//    assert(run - 1 < samples_last.size());
    auto s = sampled_tail_(run);
//    auto prev_sample = samples_last[pred_to_run[pred_pos] - 1];

    return (prev_sample + delta) % bwt_size_;
  }

 private:
  std::size_t bwt_size_;
  TPred pred_;
  TPredToRun pred_to_run_;
  TSampledTail sampled_tail_;
};

template<typename TPhi, typename TSplit, typename TLF, typename TSample>
class PhiForRange {
 public:
  template<typename TReporter, typename Range>
  std::pair<bool, std::size_t> compute(const Range &t_range,
                                       std::pair<bool, std::size_t> t_last_value,
                                       std::size_t t_n_jumps,
                                       TReporter &t_reporter) const {

    auto first = t_range.first;
    auto last = t_range.second;

    if (last < first) { return {false, -1}; }

    auto last_value = t_last_value;
    if (t_n_jumps == sampling_size_) {
      // Reach the limits of backward jumps, so the last value is valid
      last_value.first = true;
    } else {
      auto sample = sample_at_(last);
      if (sample.has_value()) {
        last_value = {true, sample.value() + t_n_jumps};
      }
    }

    // Report the last values while they are valid
    while (first <= last && last_value.first) {
      t_reporter(last_value.second);

      last_value = phi_(last_value.second);
      --last;
    }

    if (last < first) { return last_value; }

    auto runs_in_range = split_(first, last);
    for (auto it = rbegin(runs_in_range); it != rend(runs_in_range); ++it) {
      last_value = compute(lf_(*it.range, *it.c), last_value, t_n_jumps + 1, t_reporter);
    }

    return last_value;
  }

 private:
  TPhi phi_;
  TSplit split_;
  TLF lf_;
  TSample sample_at_;

  std::size_t sampling_size_;
};

}
#endif //RI_PHI_H_
