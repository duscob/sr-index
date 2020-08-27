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

template<typename TPredecessor, typename TPredecessorToTailRun, typename TSampledTail>
class Phi {
 public:
  Phi(const TPredecessor &t_predecessor,
      const TPredecessorToTailRun &t_predecessor_to_tail_run,
      const TSampledTail &t_sampled_tail,
      std::size_t t_bwt_size)
      : predecessor_{t_predecessor},
        predecessor_to_tail_run_{t_predecessor_to_tail_run},
        sampled_tail_{t_sampled_tail},
        bwt_size_{t_bwt_size} {
  }

/*
 * Phi function. Phi(SA[0]) is undefined
 */
  std::pair<std::size_t, bool> operator()(std::size_t t_prev_value) {

    assert(t_prev_value != bwt_size_ - 1);

    auto p = predecessor_(t_prev_value);

    //pred_pos is the rank of the predecessor of t_prev_value (circular)
    auto pred_pos = p.first;

    //the actual predecessor
    auto pred_value = p.second;

    //distance from predecessor
    auto delta = pred_value < t_prev_value ? t_prev_value - pred_value : t_prev_value + 1;

    //cannot fall on first head run: this can happen only if I call Phi(SA[0])
    auto run = predecessor_to_tail_run_(pred_pos);

    //sample at the end of previous run
    auto prev_sample = sampled_tail_(run.first);

    return {(prev_sample + delta) % bwt_size_, run.second};
  }

 private:
  TPredecessor predecessor_;
  TPredecessorToTailRun predecessor_to_tail_run_;
  TSampledTail sampled_tail_;

  std::size_t bwt_size_;
};

template<typename TPredecessor, typename TPredecessorToTailRun, typename TSampledTail>
auto buildPhi(const TPredecessor &t_predecessor,
              const TPredecessorToTailRun &predecessor_to_tail_run,
              const TSampledTail &t_sampled_tail,
              std::size_t t_bwt_size) {
  return Phi<TPredecessor, TPredecessorToTailRun, TSampledTail>(t_predecessor,
                                                                predecessor_to_tail_run,
                                                                t_sampled_tail,
                                                                t_bwt_size);
}

template<typename TPhi, typename TSplit, typename TBackward, typename TSample>
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
        // Position last is sampled, so we can use the sampled value
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
      last_value = compute(backward_(*it.range, *it.c), last_value, t_n_jumps + 1, t_reporter);
    }

    return last_value;
  }

 private:
  TPhi phi_;
  TSplit split_;
  TBackward backward_; // LF
  TSample sample_at_;

  std::size_t sampling_size_;
};

}
#endif //RI_PHI_H_
