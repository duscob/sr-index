//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_TOOLS_H_
#define RI_TOOLS_H_

#include <cstddef>
#include <experimental/optional>

#include "toehold.h"

namespace sri {

class GetSameArg {
 public:
  template<typename T>
  auto operator()(const T &t) const {
    return t;
  }
};

template<typename TContainer>
class RandomAccessForCRefContainer {
 public:
  explicit RandomAccessForCRefContainer(std::reference_wrapper<const TContainer> t_container)
      : container_{t_container} {}

  auto operator()(std::size_t i) const {
    return container_.get()[i];
  }

  auto operator[](std::size_t i) const {
    return (*this)(i);
  }

 private:
  std::reference_wrapper<const TContainer> container_;
};

template<typename TContainer>
class RandomAccessForContainer {
 public:
  explicit RandomAccessForContainer(const TContainer &t_container) : container_{t_container} {}

  auto operator()(std::size_t i) const {
    return container_.get()[i];
  }

 private:
  TContainer container_;
};

template<typename TContainer>
auto buildRandomAccessForContainer(const TContainer &t_container) {
  return RandomAccessForContainer<TContainer>(t_container);
}

template<typename TContainer1, typename TContainer2>
class RandomAccessForTwoContainers {
 public:
  RandomAccessForTwoContainers(const TContainer1 &t_container1, const TContainer2 &t_container2)
      : container1_{t_container1}, container2_{t_container2} {
  }

  auto operator()(std::size_t i) const {
    return std::make_pair(container1_.get()[i], container2_.get()[i]);
  }

 private:
  TContainer1 container1_;
  TContainer2 container2_;
};

template<typename TContainer1, typename TContainer2>
auto buildRandomAccessForTwoContainers(const TContainer1 &t_container1, const TContainer2 &t_container_2) {
  return RandomAccessForTwoContainers<TContainer1, TContainer2>(t_container1, t_container_2);
}

template<typename TContainer1>
class RandomAccessForTwoContainersDefault {
 public:
  RandomAccessForTwoContainersDefault(const TContainer1 &t_container1, bool t_default_value)
      : container1_{t_container1}, default_value_{t_default_value} {
  }

  auto operator()(std::size_t i) const {
    return std::make_pair(container1_.get()[i], default_value_);
  }

 private:
  TContainer1 container1_;
  bool default_value_;
};

template<typename TContainer1>
auto buildRandomAccessForTwoContainersDefault(const TContainer1 &t_container1, bool t_default_value) {
  return RandomAccessForTwoContainersDefault<TContainer1>(t_container1, t_default_value);
}

class SampleValidatorDefault {
 public:
  bool operator()(std::size_t /*t_pred_pos*/, std::size_t /*t_delta*/, bool t_is_trusted_pred) const {
    return t_is_trusted_pred;
  }
};

template<typename TGetTrustedAreaIdx, typename TGetTrustedArea>
class SampleValidator {
 public:
  SampleValidator(const TGetTrustedAreaIdx &t_get_trusted_area_idx, const TGetTrustedArea &t_get_trusted_area)
      : get_trusted_area_idx_{t_get_trusted_area_idx}, get_trusted_area_{t_get_trusted_area} {
  }

  bool operator()(std::size_t t_pred_pos, std::size_t t_delta, bool t_is_trusted_pred) const {
    if (!t_is_trusted_pred) {
      auto idx = get_trusted_area_idx_(t_pred_pos);
      auto trusted_area = get_trusted_area_(idx);
      if (trusted_area <= t_delta) {
        return false;
      }
    }

    return true;
  }

 private:
  TGetTrustedAreaIdx get_trusted_area_idx_;
  TGetTrustedArea get_trusted_area_;
};

template<typename TGetTrustedAreaIdx, typename TGetTrustedArea>
auto buildSampleValidator(const TGetTrustedAreaIdx &t_get_trusted_area_idx, const TGetTrustedArea &t_get_trusted_area) {
  return SampleValidator<TGetTrustedAreaIdx, TGetTrustedArea>(t_get_trusted_area_idx, t_get_trusted_area);
}

template<typename TRLEString>
class SplitInRuns {
 public:
  explicit SplitInRuns(const TRLEString &t_string) : string_{t_string} {}

  auto operator()(std::size_t t_first, std::size_t t_last) const {
    return string_.get().break_in_runs(std::make_pair(t_first, t_last));
  }

 private:
  TRLEString string_;
};

template<typename TRLEString>
auto buildSplitInRuns(const TRLEString &t_string) {
  return SplitInRuns<TRLEString>(t_string);
}

template<typename TRLEString>
class RankOfChar {
 public:
  explicit RankOfChar(const TRLEString &t_string) : string_{t_string} {}

  auto operator()(std::size_t t_pos, unsigned char t_char) const {
    return string_.get().rank(t_pos, t_char);
  }

 private:
  TRLEString string_;
};

template<typename TRLEString>
auto buildRankOfChar(const TRLEString &t_string) {
  return RankOfChar<TRLEString>(t_string);
}

template<typename TRLEString>
class RunOfSAPosition {
 public:
  explicit RunOfSAPosition(const TRLEString &t_string) : string_{t_string} {}

  auto operator()(std::size_t t_pos) const {
    return string_.get().get_run_of(t_pos);
  }

 private:
  TRLEString string_;
};

template<typename TRLEString>
auto buildRunOfSAPosition(const TRLEString &t_string) {
  return RunOfSAPosition<TRLEString>(t_string);
}

template<typename TRankRunTail, typename TSampledTextPosition>
class GetSampleForBWTRun {
 public:
  GetSampleForBWTRun(const TRankRunTail &t_rank_run_tail, const TSampledTextPosition &t_sampled_text_position)
      : rank_run_tail_{t_rank_run_tail}, sampled_text_position_{t_sampled_text_position} {
  }

  std::size_t operator()(std::size_t run) const {
    auto sampled_run = rank_run_tail_(run);
    return sampled_text_position_(sampled_run);
  }

 private:
  TRankRunTail rank_run_tail_;
  TSampledTextPosition sampled_text_position_;
};

template<typename TRankRunTail, typename TSampledTextPosition>
auto buildGetSampleForBWTRun(const TRankRunTail &t_rank_run_tail, const TSampledTextPosition &t_sampled_text_position) {
  return GetSampleForBWTRun<TRankRunTail, TSampledTextPosition>(t_rank_run_tail, t_sampled_text_position);
}

template<typename TRunOfSAPosition, typename TIsRunSampled, typename TSampleForBWTRun>
class GetSampleForSAPosition {
 public:
  GetSampleForSAPosition(const TRunOfSAPosition &t_run_of_sa_position,
                         const TIsRunSampled &t_is_run_sampled,
                         const TSampleForBWTRun &t_sample_for_bwt_run)
      : run_of_sa_position_{t_run_of_sa_position},
        is_run_sampled_{t_is_run_sampled},
        sample_for_bwt_run_{t_sample_for_bwt_run} {
  }

  std::experimental::optional<std::size_t> operator()(std::size_t pos) const {
    auto run = run_of_sa_position_(pos);

    if (pos != run.second || !is_run_sampled_(run.first)) { return std::experimental::nullopt; }

    return sample_for_bwt_run_(run.first);
  }

 private:
  TRunOfSAPosition run_of_sa_position_;
  TIsRunSampled is_run_sampled_;
  TSampleForBWTRun sample_for_bwt_run_;
};

template<typename TRunOfSAPosition, typename TIsRunSampled, typename TSampleForBWTRun>
auto buildGetSampleForSAPosition(const TRunOfSAPosition &t_run_of_sa_position,
                                 const TIsRunSampled &t_is_run_sampled,
                                 const TSampleForBWTRun &t_sample_for_bwt_run) {
  return GetSampleForSAPosition<TRunOfSAPosition, TIsRunSampled, TSampleForBWTRun>(
      t_run_of_sa_position, t_is_run_sampled, t_sample_for_bwt_run);
}

template<typename TRLEString>
class GetPreviousPositionInRange {
 public:
  explicit GetPreviousPositionInRange(const TRLEString &t_string) : string_{t_string} {}

  template<typename TRange, typename TChar>
  std::size_t operator()(const TRange &t_range, const TChar &t_c) const {
    //find last c in range (there must be one because range1 is not empty)
    //and get its sample (must be sampled because it is at the end of a run)
    //note: by previous check, bwt[range.second] != c, so we can use argument range.second
    auto rnk = string_.rank(t_range.second, t_c);

    //there must be at least one c before range.second
    assert(rnk > 0);

    //this is the rank of the last c
    rnk--;

    //jump to the corresponding BWT position
    auto j = string_.select(rnk, t_c);

    //the c must be in the range
    assert(t_range.first <= j && j < t_range.second);

    return j;
  }

 private:
  TRLEString string_;
};

template<typename TRLEString>
auto buildGetPreviousPositionInRange(const TRLEString &t_string) {
  return GetPreviousPositionInRange<TRLEString>(t_string);
}

template<typename TRLEString, typename TSampleForSAPosition>
class GetLastValue {
 public:
  GetLastValue(const TRLEString &t_string, const TSampleForSAPosition &t_sample_for_sa_position)
      : string_{t_string}, sample_for_sa_position_{t_sample_for_sa_position} {
  }

  template<typename TRange, typename TChar>
  std::experimental::optional<std::size_t> operator()(
      const TRange &t_range,
      const TRange &t_next_range,
      const TChar &t_c,
      std::size_t /*t_step*/,
      const std::experimental::optional<std::size_t> &t_last_value) const {

    if (t_next_range.second < t_next_range.first) { return t_last_value; }

    if (string_.get()[t_range.second] == t_c) {
      if (t_last_value) {
        // last c is at the end of range. Then, we have this sample by induction!
        assert(0 < *t_last_value);

        return *t_last_value - 1;
      } else {
        if (string_.get()[t_range.second + 1] == t_c) {
          // last position of the range is not at the end of a BWT run
          return t_last_value;
        }

        // last position of the range is also at the end of a BWT run
        return sample_for_sa_position_(t_range.second);
      }
    } else {
      //find last c in range (there must be one because range1 is not empty)
      //and get its sample (must be sampled because it is at the end of a run)
      //note: by previous check, bwt[range.second] != c, so we can use argument range.second
      auto rnk = string_.get().rank(t_range.second, t_c);

      //there must be at least one c before range.second
      assert(rnk > 0);

      //this is the rank of the last c
      rnk--;

      //jump to the corresponding BWT position
      auto j = string_.get().select(rnk, t_c);

      //the c must be in the range
      assert(t_range.first <= j && j < t_range.second);

      return sample_for_sa_position_(j);
    }
  }

 private:
  TRLEString string_;
  TSampleForSAPosition sample_for_sa_position_;
};

template<typename TRLEString, typename TSampleForSAPosition>
auto buildGetLastValue(const TRLEString &t_string, const TSampleForSAPosition &t_sample_for_sa_position) {
  return GetLastValue<TRLEString, TSampleForSAPosition>(t_string, t_sample_for_sa_position);
}

template<typename TGetSample, typename TNavigate>
class ComputeSAValue {
 public:
  ComputeSAValue(const TGetSample &t_get_sample,
                 const TNavigate &t_navigate,
                 std::size_t t_seq_size,
                 bool t_is_backward_nav)
      : get_sample_{t_get_sample}, navigate_{t_navigate}, seq_size_{t_seq_size}, is_backward_nav_{t_is_backward_nav} {
  }

  std::size_t operator()(std::size_t t_idx) const {
    std::size_t n_jumps = 0;
    auto sample = get_sample_(t_idx); // std::optional

    while (!sample) {
      t_idx = navigate_(t_idx);
      sample = get_sample_(t_idx);
      ++n_jumps;
    }

    // If the navigation is backward we go forward n_jumps jumps, else we go back n_jumps steps.
    return (*sample + 1 + seq_size_ + (is_backward_nav_ ? 1 : -1) * n_jumps) % seq_size_;
  }

 private:
  TGetSample get_sample_;
  TNavigate navigate_;

  std::size_t seq_size_;
  bool is_backward_nav_;
};

template<typename TGetSample, typename TNavigateBackward>
auto buildComputeSAValueBackward(const TGetSample &t_sample_at, const TNavigateBackward &t_lf, std::size_t t_bwt_size) {
  return ComputeSAValue<TGetSample, TNavigateBackward>(t_sample_at, t_lf, t_bwt_size, true);
}

template<typename TGetSample, typename TNavigateForward>
auto buildComputeSAValueForward(const TGetSample &t_sample_at, const TNavigateForward &t_psi, std::size_t t_bwt_size) {
  return ComputeSAValue<TGetSample, TNavigateForward>(t_sample_at, t_psi, t_bwt_size, false);
}

class GetOptionalValue {
 public:
  explicit GetOptionalValue(std::size_t t_final_value) : final_value_{t_final_value} {
  }

  auto operator()(std::size_t /*t_step*/) const {
    return std::experimental::make_optional(final_value_);
  }

 private:
  std::size_t final_value_;
};

}

#endif //RI_TOOLS_H_
