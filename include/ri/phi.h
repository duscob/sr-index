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

template<typename TPredecessor, typename TPredecessorToTailRun, typename TSampledTail, typename TSampledTailValidator>
class Phi {
 public:
  Phi(const TPredecessor &t_predecessor,
      const TPredecessorToTailRun &t_predecessor_to_tail_run,
      const TSampledTail &t_sampled_tail,
      const TSampledTailValidator &t_is_sampled_tail_valid,
      std::size_t t_bwt_size)
      : predecessor_{t_predecessor},
        predecessor_to_tail_run_{t_predecessor_to_tail_run},
        sampled_tail_{t_sampled_tail},
        is_sampled_tail_valid_{t_is_sampled_tail_valid},
        bwt_size_{t_bwt_size} {
  }

/*
 * Phi function. Phi(SA[0]) is undefined
 */
  std::pair<std::size_t, bool> operator()(std::size_t t_prev_value) const {

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

    auto sampled_tail_validity = is_sampled_tail_valid_(pred_pos, delta, run.second);
    return {(prev_sample + delta) % bwt_size_, sampled_tail_validity};
  }

 private:
  TPredecessor predecessor_;
  TPredecessorToTailRun predecessor_to_tail_run_;
  TSampledTail sampled_tail_;
  TSampledTailValidator is_sampled_tail_valid_;

  std::size_t bwt_size_;
};

template<typename TPredecessor, typename TPredecessorToTailRun, typename TSampledTail, typename TSampledTailValidator>
auto buildPhi(const TPredecessor &t_predecessor,
              const TPredecessorToTailRun &predecessor_to_tail_run,
              const TSampledTail &t_sampled_tail,
              const TSampledTailValidator &t_sampled_tail_validator,
              std::size_t t_bwt_size) {
  return Phi<TPredecessor, TPredecessorToTailRun, TSampledTail, TSampledTailValidator>(
      t_predecessor, predecessor_to_tail_run, t_sampled_tail, t_sampled_tail_validator, t_bwt_size);
}

template<typename TPhi, typename TSplitInBWTRun, typename TBackwardNav, typename TSampleAt>
class PhiForRange {
 public:
  PhiForRange(const TPhi &t_phi,
              const TSplitInBWTRun &t_split,
              const TBackwardNav &t_lf,
              const TSampleAt &t_sample_at,
              std::size_t t_sampling_size,
              std::size_t t_bwt_size)
      : phi_{t_phi},
        split_{t_split},
        lf_{t_lf},
        sample_at_{t_sample_at},
        sampling_size_{t_sampling_size},
        bwt_size_{t_bwt_size} {
  }

  template<typename TReport, typename TRange>
  void operator()(const TRange &t_range, std::size_t t_prev_value, TReport &t_report) const {
    auto last_value = phi_(t_prev_value);
    compute(t_range, last_value, 0, t_report);
  }

  template<typename TReporter, typename Range>
  std::pair<std::size_t, bool> compute(const Range &t_range,
                                       std::pair<std::size_t, bool> t_last_value,
                                       std::size_t t_n_jumps,
                                       TReporter &t_reporter) const {

    auto first = t_range.first;
    auto last = t_range.second;

    if (last < first) { return {false, -1}; }

    auto last_value = t_last_value;
    if (sampling_size_ <= t_n_jumps) {
      // Reach the limits of backward jumps, so the last value is valid
      last_value.second = true;
    }

    do {

      // Report the last values while they are valid
      while (first <= last && last_value.second) {
        t_reporter(last_value.first);

        last_value = phi_(last_value.first);
        --last;
      }

      if (last < first) { return last_value; }

      auto sample = sample_at_(last);
      if (sample) {
        // Position last is sampled, so we can use the sampled value

        // NOTE we sampled the position of the i-th BWT char, but here we want the position of SA[i], so we need + 1
        last_value = {(sample.value() + 1 + t_n_jumps) % bwt_size_, true};
      }

    } while (last_value.second);

    // TODO Don't split it in all sub-runs. We can go in depth recursively with only the last remaining sub-run, and continue with the other at the same level.
    auto runs_in_range = split_(first, last);

    auto it = rbegin(runs_in_range);
    last_value = compute(lf_(it->range, it->c), last_value, t_n_jumps + 1, t_reporter);
    while(++it != rend(runs_in_range)) {
      last_value = compute(it->range, last_value, t_n_jumps, t_reporter);
    }

    return last_value;
  }

 private:
  TPhi phi_;
  TSplitInBWTRun split_; // Split an interval in its internal BWT runs
  TBackwardNav lf_; // LF
  TSampleAt sample_at_; // Access to last value of a BWT run. Note that some tails are not sampled.

  std::size_t sampling_size_;
  std::size_t bwt_size_;
};

template<typename TPhi, typename TSplitInBWTRun, typename TBackwardNav, typename TSampleAt>
auto buildPhiForRange(const TPhi &t_phi,
                      const TSplitInBWTRun &t_split,
                      const TBackwardNav &t_lf,
                      const TSampleAt &t_sample_at,
                      std::size_t t_sampling_size,
                      std::size_t t_bwt_size) {
  return PhiForRange<TPhi, TSplitInBWTRun, TBackwardNav, TSampleAt>(
      t_phi, t_split, t_lf, t_sample_at, t_sampling_size, t_bwt_size);
}

template<typename TPhi>
class ComputeAllValuesWithPhi {
 public:
  explicit ComputeAllValuesWithPhi(const TPhi &t_phi) : phi_{t_phi} {}

  template<typename TRange, typename TReport>
  void operator()(const TRange &t_range, std::experimental::optional<std::size_t> t_k, TReport &t_report) const {
    auto k = *t_k;
    for (auto i = t_range.first; i <= t_range.second; ++i) {
      t_report(k);

      // This function is for original r-index, so phi returns always a valid value
      k = phi_(k).first;
    }
  }

 private:
  TPhi phi_;
};

template<typename TPhi>
auto buildComputeAllValuesWithPhi(const TPhi &t_phi) {
  return ComputeAllValuesWithPhi<TPhi>(t_phi);
}

template<typename TPhiForRange, typename TGetValueForSAPosition>
class ComputeAllValuesWithPhiForRange {
 public:
  ComputeAllValuesWithPhiForRange(const TPhiForRange &t_phi_for_range,
                                  const TGetValueForSAPosition &t_get_value_for_sa_position)
      : phi_for_range_{t_phi_for_range}, get_value_for_sa_position_{t_get_value_for_sa_position} {
  }

  template<typename TRange, typename TReport>
  void operator()(const TRange &t_range, std::experimental::optional<std::size_t> t_k, TReport &t_report) const {
    auto k = t_k ? *t_k : get_value_for_sa_position_(t_range.second);
    t_report(k);

    auto range = t_range;
    --range.second;

    phi_for_range_(range, k, t_report);
  }

 private:
  TPhiForRange phi_for_range_;
  TGetValueForSAPosition get_value_for_sa_position_;
};

template<typename TPhiForRange, typename TGetValueForSAPosition>
auto buildComputeAllValuesWithPhiForRange(const TPhiForRange &t_phi_for_range,
                                          const TGetValueForSAPosition &t_get_value_for_sa_position) {
  return ComputeAllValuesWithPhiForRange<TPhiForRange, TGetValueForSAPosition>(
      t_phi_for_range, t_get_value_for_sa_position);
}

}
#endif //RI_PHI_H_
