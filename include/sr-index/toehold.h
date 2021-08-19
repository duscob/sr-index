//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/17/21.
//

#ifndef SRI_TOEHOLD_H_
#define SRI_TOEHOLD_H_

#include <functional>

#include "psi.h"

namespace sri {

//! Data for the current backward-search step to compute toehold value for final range
//! \tparam TChar
template<typename TChar>
struct DataBackwardSearchStep {
  TChar c; // Character for LF step in the range
  std::size_t step; // Number of LF steps to reach final range
  std::size_t range_start; // Start of range before LF step
  std::size_t range_end; // End of range before LF step
};

//! Compute data corresponding to backward-search step for given range and character
//! \tparam TIsLFTrivial Predicate to check if the LF step in the range with the character is trivial, i.e. if lf(range, c) matches with next_range in the corresponding limit.
template<typename TIsLFTrivial>
class GetLastSpecialBackwardSearchStep {
 public:
  explicit GetLastSpecialBackwardSearchStep(const TIsLFTrivial &t_is_lf_trivial) : is_lf_trivial_{t_is_lf_trivial} {
  }

  template<typename TRange, typename TChar, typename TDataBackwardSearchStep>
  auto operator()(
      const TRange &t_range,
      const TRange &t_next_range,
      const TChar &t_c,
      std::size_t t_step,
      const TDataBackwardSearchStep &t_last_special_step) const {
    // TODO It is required t_next_range? If next range is empty the pattern was not found.
    // TODO Replace TDataBackwardSearchStep by DataBackwardSearchStep<TChar>

    if (t_next_range.second < t_next_range.first || is_lf_trivial_(t_range, t_c)) {
      // Next range is empty or
      // LF step on a current range limit (range end for backward phi or range start for forward phi) goes to next range limit.
      return t_last_special_step;
    }

    return TDataBackwardSearchStep{t_c, t_step, t_range.first, t_range.second};
  }

 private:
  TIsLFTrivial is_lf_trivial_;
};

template<typename TRLEString>
auto buildGetLastSpecialBackwardSearchStepForPhiBackward(const std::reference_wrapper<const TRLEString> &t_bwt) {
  auto is_lf_trivial_with_bwt = [t_bwt](const auto &tt_range, const auto &tt_c) {
    return t_bwt.get()[tt_range.second] == tt_c;
  };

  return GetLastSpecialBackwardSearchStep(is_lf_trivial_with_bwt);
}

template<typename TBitVector>
auto buildGetLastSpecialBackwardSearchStepForPhiForward(const std::reference_wrapper<const PsiCore<TBitVector>> &t_psi) {
  auto is_lf_trivial_with_psi = [t_psi](const auto &tt_range, const auto &tt_c) {
    return t_psi.get().partial_psi[tt_c][tt_range.first] == 1;
  };

  return GetLastSpecialBackwardSearchStep(is_lf_trivial_with_psi);
}

//! Compute toehold value for PhiBackward using BWT rank and select.
template<typename TBWTRank, typename TBWTSelect, typename TGetSAValue>
class ComputeToeholdValueForPhiBackward {
 public:
  ComputeToeholdValueForPhiBackward(
      std::size_t t_n, const TBWTRank &t_bwt_rank, const TBWTSelect &t_bwt_select, const TGetSAValue &t_get_sa_value)
      : n_{t_n}, bwt_rank_{t_bwt_rank}, bwt_select_{t_bwt_select}, get_sa_value_{t_get_sa_value} {
  }

  template<typename TChar, typename TRange>
  auto operator()(const DataBackwardSearchStep<TChar> &t_data, const TRange &) const {
    // Find last c in range [t_data.range_start..t_data.range_end] (there must be one because final range is not empty)
    // and get its sample (must be sampled because it is at the end of a run).

    // Notice that it could bwt[t_data.range_end] == c, so we must use t_data.range_end + 1 as argument of bwt rank.
    auto rnk = bwt_rank_(t_data.range_end + 1, t_data.c);

    //there must be at least one c before range.second
    assert(rnk > 0);

    //this is the rank of the last c
    --rnk;

    //jump to the corresponding BWT position
    auto j = bwt_select_(rnk, t_data.c);

    //the c must be in the range
    assert(t_data.range_start <= j && j <= t_data.range_end);

    return (get_sa_value_(j) + n_ - t_data.step - 1) % n_;
  }

 private:
  std::size_t n_;

  TBWTRank bwt_rank_;
  TBWTSelect bwt_select_;

  TGetSAValue get_sa_value_;
};

template<typename TRLEString, typename TGetSAValue>
auto buildComputeToeholdValueForPhiBackward(const std::reference_wrapper<const TRLEString> &t_bwt,
                                            const TGetSAValue &t_get_sa_value) {
  auto bwt_rank = [t_bwt](auto tt_pos, auto tt_c) { return t_bwt.get().rank(tt_pos, tt_c); };
  auto bwt_select = [t_bwt](auto tt_rnk, auto tt_c) { return t_bwt.get().select(tt_rnk, tt_c); };

  return ComputeToeholdValueForPhiBackward(t_bwt.get().size(), bwt_rank, bwt_select, t_get_sa_value);
}

}

#endif //SRI_TOEHOLD_H_
