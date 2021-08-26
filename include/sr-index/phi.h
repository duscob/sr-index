//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef SRI_PHI_H_
#define SRI_PHI_H_

#include <cstddef>
#include <cassert>
#include <iterator>
#include <experimental/optional>

namespace sri {

//! Get the sample associated to a marked value position.
/**
 * @tparam TGetSampleIndex Get the index of the sample associated to the given marked value.
 * @tparam TGetSample Get the sample value for the given index.
 * @tparam TSampleValidator Validate the given sample.
 */
template<typename TGetSampleIndex, typename TGetSample, typename TSampleValidator>
class GetSampleForPhi {
 public:
  GetSampleForPhi(const TGetSampleIndex &t_get_sample_index,
                  const TGetSample &t_get_sample,
                  const TSampleValidator &t_is_sample_valid)
      : get_sample_index_{t_get_sample_index}, get_sample_{t_get_sample}, is_sample_valid_{t_is_sample_valid} {
  }

  //! Get sample
  /**
   * @param t_marked_value_pos Position of the marked value.
   * @param t_delta Difference between original value and marked value.
   * @return <s, v>: sample value and its validity.
   */
  auto operator()(std::size_t t_marked_value_pos, std::size_t t_delta) const {
    //cannot fall on first head run: this can happen only if I call PhiBackward(SA[0])
    auto[sample_index, initial_sample_validity] = get_sample_index_(t_marked_value_pos);

    auto sample = get_sample_(sample_index);

    auto sample_validity = is_sample_valid_(t_marked_value_pos, t_delta, initial_sample_validity);

    return std::make_pair(sample, sample_validity);
  }

 private:
  TGetSampleIndex get_sample_index_;
  TGetSample get_sample_;
  TSampleValidator is_sample_valid_;
};

//! Phi function
/**
 * @param t_prev_value Previous value of the Suffix Array.
 * @param t_get_related_value Get value related to given SA value.
 * @param t_get_sample Get sample associated to related value.
 * @param t_get_delta Get difference between SA given value and related marked value.
 * @param t_get_next_value Get next value using sample and delta.
 * @return <s, v>: next value in SA and its validity.
 */
template<typename TGetRelatedValue, typename TGetSample, typename TGetDelta, typename TGetNextValue>
auto phi(std::size_t t_prev_value,
         const TGetRelatedValue &t_get_related_value,
         const TGetSample &t_get_sample,
         const TGetDelta &t_get_delta,
         const TGetNextValue &t_get_next_value) {
  // { Related value position/rank (circular); Actual related value }
  auto[related_value_pos, related_value] = t_get_related_value(t_prev_value);

  // Distance between the previous value and its related value
  auto delta = t_get_delta(t_prev_value, related_value);

  // Sample associated with related value and its validity
  auto[sample, sample_validity] = t_get_sample(related_value_pos, delta);

  return std::make_pair(t_get_next_value(sample, delta), sample_validity);
}

//! Phi function based on backward search using predecessor data structure.
/**
 * @tparam TPredecessor Get predecessor (circular) for a given value.
 * @tparam TGetSample Get sample associated to related value found.
 */
template<typename TPredecessor, typename TGetSample>
class PhiBackward {
 public:
  PhiBackward(const TPredecessor &t_predecessor, const TGetSample &t_get_sample, std::size_t t_size)
      : predecessor_{t_predecessor}, get_sample_{t_get_sample}, size_{t_size} {
  }

  //! Phi backward function. PhiBackward(SA[0]) is undefined.
  /**
   * @param t_prev_value Previous value of the Suffix Array.
   * @return <s, v>: next value in SA and its validity.
   * @note Cannot fall on first head run: this can happen only if I call PhiBackward(SA[0])
   */
  auto operator()(std::size_t t_prev_value) const {
    assert(t_prev_value != this->size_ - 1);

    // Distance between the previous value and its related value
    auto get_delta = [](std::size_t t_prev_value, std::size_t t_related_value) {
      return t_related_value < t_prev_value ? t_prev_value - t_related_value : t_prev_value + 1;
    };

    // Next value using sample and delta
    auto get_next_value = [size = this->size_](std::size_t t_sample, std::size_t t_delta) {
      return (t_sample + t_delta) % size;
    };

    return phi(t_prev_value, predecessor_, get_sample_, get_delta, get_next_value);
  }

 private:
  TPredecessor predecessor_; // Smaller than...
  TGetSample get_sample_;

  std::size_t size_;
};

template<typename TPredecessor, typename TGetSampleIndex, typename TGetSample, typename TSampleValidator>
auto buildPhiBackward(const TPredecessor &t_predecessor,
                      const TGetSampleIndex &t_get_sample_index,
                      const TGetSample &t_get_sample,
                      const TSampleValidator &t_sample_validator,
                      std::size_t t_bwt_size) {
  auto get_sample = GetSampleForPhi(t_get_sample_index, t_get_sample, t_sample_validator);
  return PhiBackward<TPredecessor, decltype(get_sample)>(t_predecessor, get_sample, t_bwt_size);
}

//! Phi function based on forward search using soft-successor data structure.
/**
 * @tparam TSoftSuccessor Get soft-successor (circular) for a given value. The returned value is greater or equal than the given value.
 * @tparam TGetSample Get sample associated to related value found.
 */
template<typename TSoftSuccessor, typename TGetSample>
class PhiForward {
 public:
  PhiForward(const TSoftSuccessor &t_soft_successor, const TGetSample &t_get_sample, std::size_t t_size)
      : soft_successor_{t_soft_successor}, get_sample_{t_get_sample}, size_{t_size} {
  }

  //! Phi forward function. PhiForward(SA[n-1]) is undefined.
  /**
   * @param t_prev_value Previous value of the Suffix Array.
   * @return <s, v>: next value in SA and its validity.
   * @note Cannot fall on last tail run: this can happen only if I call PhiForward(SA[n-1])
   */
  std::pair<std::size_t, bool> operator()(std::size_t t_prev_value) const {
    // Distance between the previous value and its related value
    auto get_delta = [](std::size_t t_prev_value, std::size_t t_related_value) {
      return t_related_value - t_prev_value;
    };

    // Next value using sample and delta
    auto get_next_value = [](std::size_t t_sample, std::size_t t_delta) {
      return t_sample - t_delta;
    };

    return phi(t_prev_value, soft_successor_, get_sample_, get_delta, get_next_value);
  }

 private:
  TSoftSuccessor soft_successor_; // Greater or equal than...
  TGetSample get_sample_;

  std::size_t size_;
};

template<typename TSoftSuccessor, typename TGetSampleIndex, typename TGetSample, typename TSampleValidator>
auto buildPhiForward(const TSoftSuccessor &t_soft_successor,
                     const TGetSampleIndex &t_get_sample_index,
                     const TGetSample &t_get_sample,
                     const TSampleValidator &t_sample_validator,
                     std::size_t t_bwt_size) {
  auto get_sample = GetSampleForPhi(t_get_sample_index, t_get_sample, t_sample_validator);
  return PhiForward<TSoftSuccessor, decltype(get_sample)>(t_soft_successor, get_sample, t_bwt_size);
}

//! Compute the links between marked positions (BWT tails) and its corresponding samples (BWT heads) to be used in PhiForward function.
/**
 * @param t_n Number of symbols in the sequence
 * @param t_r Number of BWT runs
 * @param t_select_mark_sa_pos Select function over position of marked values (BWT run tails)
 * @param t_rank_sample_sa_pos Rank function over position of sampled values (BWT run heads)
 * @param t_lf LF function
 * @param t_psi Psi function
 * @param t_report Reporter for the links
 */
template<typename TSelectMarkSAPos, typename TRankSampleSAPos, typename TLF, typename TPsi, typename TReport>
void computeMarkToSampleLinkForPhiForward(std::size_t t_n,
                                          std::size_t t_r,
                                          const TSelectMarkSAPos &t_select_mark_sa_pos,
                                          const TRankSampleSAPos &t_rank_sample_sa_pos,
                                          const TLF &t_lf,
                                          const TPsi &t_psi,
                                          TReport &t_report) {
  auto rank_sample = [&t_rank_sample_sa_pos](const auto &tt_k) {
    return t_rank_sample_sa_pos(tt_k + 1) - 1;
  };

  for (std::size_t i = 1; i <= t_r; ++i) {
    t_report(i - 1, computeMarkToSampleLinkForPhiForward(t_select_mark_sa_pos(i), t_n, t_lf, t_psi, rank_sample));
  }
}

//! Compute link between a marked position (BWT tail) and its corresponding sample (BWT heads) to be used in PhiForward function.
/**
 * @param t_mark_pos Marked position or position of the last letter of a BWT run
 * @param t_n Number of symbols in the sequence
 * @param t_lf LF function
 * @param t_psi Psi function
 * @param t_rank_sample Compute index (rank) for sample position in BWT
 * @param t_report Reporter for the links
 * @return Index of sample associated to given mark position
 */
template<typename TLF, typename TPsi, typename TRankSample>
auto computeMarkToSampleLinkForPhiForward(std::size_t t_mark_pos,
                                          std::size_t t_n,
                                          const TLF &t_lf,
                                          const TPsi &t_psi,
                                          const TRankSample &t_rank_sample) {
  // Current BWT run tail could travel together its following symbol until previous LF step
  auto j = t_lf(t_mark_pos); // Position of previous symbol

  // Psi value of the following symbol (ISA[SA[j + 1] + 1]) is the BWT run head (sample) associated to current BWT run tail (mark value)
  auto k = t_psi((j + 1) % t_n); // Psi position of next symbol in SA

  return t_rank_sample(k);
}

template<typename TPhi, typename TSplitInBWTRun, typename TBackwardNav, typename TSampleAt>
class PhiForRangeSimple {
 public:
  PhiForRangeSimple(const TPhi &t_phi,
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
//    auto last_value = phi_(t_prev_value);
    auto last_value = std::make_pair(t_prev_value, false);
    compute(t_range, last_value, 0, t_report);
  }

  template<typename TReporter, typename Range>
  std::pair<std::size_t, bool> compute(const Range &t_range,
                                       std::pair<std::size_t, bool> t_last_value,
                                       std::size_t t_n_jumps,
                                       TReporter &t_reporter) const {

    auto first = t_range.first;
    auto last = t_range.second;

    if (last < first) { return {-1, false}; }

    auto last_value = t_last_value;
    if (sampling_size_ <= t_n_jumps) {
      // Reach the limits of backward jumps, so phi for the previous value is valid
      last_value = phi_(last_value.first);

      do {
        t_reporter(last_value.first);

        last_value = phi_(last_value.first);
        --last;
      } while (first <= last);

      return last_value;
    }

    auto sample = sample_at_(last);
    if (sample) {
      // Position last is sampled, so we can use the sampled value

      // NOTE we sampled the position of the i-th BWT char, but here we want the position of SA[i], so we need + 1
      last_value = {(sample.value() + 1 + t_n_jumps) % bwt_size_, false};
      t_reporter(last_value.first);
      --last;
    }

    if (last < first) { return last_value; }

    // TODO Don't split it in all sub-runs. We can go in depth recursively with only the last remaining sub-run, and continue with the other at the same level.
    auto runs_in_range = split_(first, last);

    auto it = rbegin(runs_in_range);
    last_value = compute(lf_(it->range, it->c), last_value, t_n_jumps + 1, t_reporter);
    while (++it != rend(runs_in_range)) {
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
auto buildPhiForRangeSimple(const TPhi &t_phi,
                            const TSplitInBWTRun &t_split,
                            const TBackwardNav &t_lf,
                            const TSampleAt &t_sample_at,
                            std::size_t t_sampling_size,
                            std::size_t t_bwt_size) {
  return PhiForRangeSimple<TPhi, TSplitInBWTRun, TBackwardNav, TSampleAt>(
      t_phi, t_split, t_lf, t_sample_at, t_sampling_size, t_bwt_size);
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

    if (last < first) { return {-1, false}; }

    auto last_value = t_last_value;
    if (sampling_size_ <= t_n_jumps) {
      // Reach the limits of backward jumps, so the last value is valid
      do {
        t_reporter(last_value.first);

        last_value = phi_(last_value.first);
        --last;
      } while (first <= last);

      return last_value;
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
    while (++it != rend(runs_in_range)) {
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

template<typename TPhiForRange, typename TComputeToehold, typename TUpdateRange>
class ComputeAllValuesWithPhiForRange {
 public:
  ComputeAllValuesWithPhiForRange(const TPhiForRange &t_phi_for_range,
                                  const TComputeToehold &t_compute_toehold,
                                  const TUpdateRange &t_update_range)
      : phi_for_range_{t_phi_for_range}, compute_toehold_{t_compute_toehold}, update_range_{t_update_range} {
  }

  template<typename TRange, typename TDataLastValue, typename TReport>
  void operator()(const TRange &t_range, const TDataLastValue &t_k, TReport &t_report) const {
    // TODO In the case we need to backward search the value for the last position in the range,
    //  we can take advantage of this travel for the other positions in the same BWT sub-run
    auto k = compute_toehold_(t_k);
    t_report(k);

    auto range = update_range_(t_range);

    phi_for_range_(range, k, t_report);
  }

 private:
  TPhiForRange phi_for_range_;
  TComputeToehold compute_toehold_;
  TUpdateRange update_range_;
};

template<typename TPhiBackwardForRange, typename TGetSAValue>
auto buildComputeAllValuesWithPhiBackwardForRange(const TPhiBackwardForRange &t_phi_backward_for_range,
                                                  const TGetSAValue &t_get_sa_value) {
  auto get_new_range = [](auto tt_range) {
    --tt_range.second;
    return tt_range;
  };

  return ComputeAllValuesWithPhiForRange(t_phi_backward_for_range, t_get_sa_value, get_new_range);
}

template<typename TPhiForwardForRange, typename TGetSAValue>
auto buildComputeAllValuesWithPhiForwardForRange(const TPhiForwardForRange &t_phi_forward_for_range,
                                                 const TGetSAValue &t_get_sa_value) {
  auto get_new_range = [](auto tt_range) {
    ++tt_range.first;
    return tt_range;
  };

  return ComputeAllValuesWithPhiForRange(t_phi_forward_for_range, t_get_sa_value, get_new_range);
}

}
#endif //SRI_PHI_H_
