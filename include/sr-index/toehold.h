//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/17/21.
//

#ifndef SRI_TOEHOLD_H_
#define SRI_TOEHOLD_H_

namespace sri {

template<typename TChar>
struct DataBackwardSearchStep {
  TChar c;
  std::size_t step;
  std::size_t range_start;
  std::size_t range_end;
};

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
auto buildGetLastSpecialBackwardSearchStepForPhiBackward(const std::reference_wrapper<TRLEString> &t_bwt) {
  auto is_lf_trivial_with_bwt = [t_bwt](const auto &tt_range, const auto &tt_c) {
    return t_bwt.get()[tt_range.second] == tt_c;
  };

  return GetLastSpecialBackwardSearchStep(is_lf_trivial_with_bwt);
}

}

#endif //SRI_TOEHOLD_H_
