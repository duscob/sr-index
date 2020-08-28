//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_TOOLS_H_
#define RI_TOOLS_H_

#include <cstddef>
#include <experimental/optional>

namespace ri {

template<typename TContainer>
class RandomAccessForContainer {
 public:
  RandomAccessForContainer(const TContainer &t_container) : container_{t_container} {}

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

template<typename TRLEString>
class SplitInRuns {
 public:
  SplitInRuns(const TRLEString &t_string) : string_{t_string} {}

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
  RankOfChar(const TRLEString &t_string) : string_{t_string} {}

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
  RunOfSAPosition(const TRLEString &t_string) : string_{t_string} {}

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

template<typename TRLEString, typename TIsRunSampled, typename TSampleForBWTRun>
auto buildGetSampleForSAPosition(const TRLEString &t_string,
                                 const TIsRunSampled &t_is_run_sampled,
                                 const TSampleForBWTRun &t_sample_for_bwt_run) {
  return GetSampleForSAPosition<TRLEString, TIsRunSampled, TSampleForBWTRun>(
      t_string, t_is_run_sampled, t_sample_for_bwt_run);
}

}

#endif //RI_TOOLS_H_
