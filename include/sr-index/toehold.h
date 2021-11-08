//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/17/21.
//

#ifndef SRI_TOEHOLD_H_
#define SRI_TOEHOLD_H_

#include <cstddef>
#include <functional>

namespace sri {

//! Data for the current backward-search step to compute toehold value for final range
template<typename TRunData>
struct DataBackwardSearchStep {
  std::size_t step; // Number of LF steps to reach final range
  TRunData run_data; // Data required to compute toehold

  DataBackwardSearchStep(std::size_t t_step, const TRunData &t_run_data) : step{t_step}, run_data{t_run_data} {}
};

template<typename TChar>
struct RunDataBackward {
  TChar c; // Character for LF step in the range
  std::size_t end; // End of range before LF step

  RunDataBackward(TChar t_c, std::size_t t_end) : c{t_c}, end{t_end} {}
};

template<typename TChar>
using DataBackwardSearchStepBackward = DataBackwardSearchStep<RunDataBackward<TChar>>;

template<typename TChar>
class GetInitialDataBackwardSearchStep {
 public:
  GetInitialDataBackwardSearchStep(const TChar &t_c, std::size_t t_range_end) : c{t_c}, range_end{t_range_end} {}

  auto operator()(std::size_t t_step) const {
    return DataBackwardSearchStepBackward<TChar>{t_step, {c, range_end}};
  }

 private:
  TChar c;
  std::size_t range_end;
};

template<typename TChar>
auto buildGetInitialDataBackwardSearchStep(const TChar &t_c, std::size_t t_range_end) {
  return GetInitialDataBackwardSearchStep<TChar>(t_c, t_range_end);
}

using DataBackwardSearchStepForward = DataBackwardSearchStep<std::size_t>;

//! Compute data corresponding to backward-search step for given range and character
//! \tparam TIsLFTrivial Predicate to check if the LF step in the range with the character is trivial,
//!     i.e., if lf(range, c) matches with next_range in the corresponding limit or if the next range is empty.
template<typename TIsLFTrivial, typename TCreateData>
class ComputeDataBackwardSearchStep {
 public:
  ComputeDataBackwardSearchStep(const TIsLFTrivial &t_is_lf_trivial, const TCreateData &t_create_data)
      : is_lf_trivial_{t_is_lf_trivial}, create_data_{t_create_data} {}

  template<typename TRange, typename TNextRange, typename TChar, typename TRunData>
  auto operator()(const TRange &t_range,
                  const TNextRange &t_next_range,
                  const TChar &t_c,
                  std::size_t t_step,
                  const DataBackwardSearchStep<TRunData> &t_last_special_step) const {
    if (is_lf_trivial_(t_range, t_c, t_next_range)) {
      // Next range is empty or
      // LF step on a current range limit (range end for backward phi or range start for forward phi) goes to next range limit.
      return t_last_special_step;
    }

    return create_data_(t_range, t_c, t_next_range, t_step);
  }

 private:
  TIsLFTrivial is_lf_trivial_;
  TCreateData create_data_;
};

template<typename TRLEString>
auto buildComputeDataBackwardSearchStepForPhiBackward(const std::reference_wrapper<const TRLEString> &t_bwt) {
  auto is_lf_trivial_with_bwt = [t_bwt](const auto &tt_range, const auto &tt_c, const auto &tt_next_range) {
    const auto &[next_start, next_end] = tt_next_range;
    return next_end < next_start || t_bwt.get()[tt_range.second] == tt_c;
  };

  auto create_data = [](const auto &tt_range, const auto &tt_c, const auto &tt_next_range, const auto &tt_step) {
    const auto &[start, end] = tt_range;
    return DataBackwardSearchStep{tt_step, RunDataBackward{tt_c, end}};
  };

  return ComputeDataBackwardSearchStep(is_lf_trivial_with_bwt, create_data);
}

template<typename TCreateData>
auto buildComputeDataBackwardSearchStepForPhiForward(const TCreateData &t_create_data) {
  auto is_lf_trivial_with_psi = [](const auto &tt_range, const auto &tt_c, const auto &tt_next_range) {
    const auto &[next_start, next_end] = tt_next_range;
    const auto &[start, end] = tt_range;
    return !(next_start < next_end) || (!(start < next_start.run.start) && start < next_start.run.end);
  };

  return ComputeDataBackwardSearchStep(is_lf_trivial_with_psi, t_create_data);
}

//TODO Remove this function if is not used
auto buildComputeDataBackwardSearchStepForPhiForward() {
  auto create_data = [](const auto &tt_range, const auto &tt_c, const auto &tt_next_range, const auto &tt_step) {
    const auto &[start, end] = tt_next_range;
    return DataBackwardSearchStepForward{tt_step, start.run.start};
  };

  return buildComputeDataBackwardSearchStepForPhiForward(create_data);
}

//! Compute toehold value for Phi Backward/Forward using BWT rank and select.
template<typename TGetSAValue>
class ComputeToehold {
 public:
  ComputeToehold(const TGetSAValue &t_get_sa_value, std::size_t t_n) : get_sa_value_{t_get_sa_value}, n_{t_n} {}

  //! Find first/last symbol c in range (there must be one because final range is not empty)
  //! and get its sample (there must be sampled because it is at the start/end of a run).
  //! \param t_data Data associated to last special (no trivial) backward search step.
  //! \return Toehold value for final range.
  template<typename TRunData>
  auto operator()(const DataBackwardSearchStep<TRunData> &t_data) const {
    const auto &[n_steps, run_data] = t_data;
    return (get_sa_value_(run_data) + n_ - n_steps - 1) % n_;
  }

 private:
  TGetSAValue get_sa_value_;
  std::size_t n_;
};

template<typename TRLEString, typename TGetSAValue>
auto buildComputeToeholdForPhiBackward(const std::reference_wrapper<const TRLEString> &t_bwt,
                                       const TGetSAValue &t_get_sa_value) {
  auto compute_sa_value = [t_bwt, t_get_sa_value](const auto &tt_run_data) {
    const auto &[c, run_end] = tt_run_data;

    // Note that bwt[ep] == c could be true, so we must use ep + 1 as argument of bwt rank.
    auto rnk = t_bwt.get().rank(run_end + 1, c); // Rank of the last symbol c in range
    assert(("There must be at least one symbol c in the range", rnk > 0));

    // Note that for rle_string.select, rnk starts in 0.
    auto j = t_bwt.get().select(rnk - 1, c); // Corresponding BWT position for last symbol c in range
    assert(("Symbol c must be in the range", j <= run_end));

    return t_get_sa_value(j);
  };

  return ComputeToehold(compute_sa_value, t_bwt.get().size());
}

template<typename TGetSAValue>
auto buildComputeToeholdForPhiForward(const TGetSAValue &t_get_sa_value, std::size_t t_n) {
  return ComputeToehold(t_get_sa_value, t_n);
}

}

#endif //SRI_TOEHOLD_H_
