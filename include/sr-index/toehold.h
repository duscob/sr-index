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
template<typename TRLEString, typename TGetValueForSAPosition>
class ComputeToeholdValueForPhiBackward {
 public:
  ComputeToeholdValueForPhiBackward(const TRLEString &t_bwt,
                                    const TGetValueForSAPosition &t_get_value_for_sa_position)
      : bwt_{t_bwt}, get_value_for_sa_position_{t_get_value_for_sa_position} {
  }

  template<typename TChar, typename TRange>
  auto operator()(const DataBackwardSearchStep<TChar> &t_data, const TRange &t_range) const {
    //find last c in range (there must be one because range1 is not empty)
    //and get its sample (must be sampled because it is at the end of a run)
    //note: by previous check, bwt[range.second] != c, so we can use argument range.second
    auto rnk = bwt_.get().rank(t_data.range_end, t_data.c); // Fixme It is posible that bwt[t_data.range_second] == t_data.c, so previous condition could be false and we need bwt_rank(t_data.range_end + 1, t_data.c)

    //there must be at least one c before range.second
    assert(rnk > 0);

    //this is the rank of the last c
    rnk--;

    //jump to the corresponding BWT position
    auto j = bwt_.get().select(rnk, t_data.c);

    //the c must be in the range
    assert(t_range.first <= j && j < t_range.second); // FIXME False since t_data.step >= 0, so we need LF steps t_data.step-times to reach t_range

    return (get_value_for_sa_position_(j) + bwt_.get().size() - t_data.step - 1) % bwt_.get().size();
  }

 private:
  TRLEString bwt_;
  TGetValueForSAPosition get_value_for_sa_position_;
};

template<typename TRLEString, typename TGetValueForSAPosition>
auto buildComputeFinalValueWithLastSpecialBackwardSearchStep(
    const TRLEString &t_bwt, const TGetValueForSAPosition &t_get_value_for_sa_position) {
  return ComputeToeholdValueForPhiBackward<TRLEString, TGetValueForSAPosition>(
      t_bwt, t_get_value_for_sa_position);
}

}

#endif //SRI_TOEHOLD_H_
