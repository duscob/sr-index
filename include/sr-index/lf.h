//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/20/21.
//

#ifndef SRI_LF_H_
#define SRI_LF_H_

namespace sri {

//! LF function
//! \tparam TRankC Rank function for character @p c
//! \tparam TCumulativeC Function to cumulative count for alphabet
//! \tparam TValue Range position type
//! \tparam TChar Char type
//! \tparam TCreateNewRange Function to create new range
//! \tparam TNewRange New range type
//! \param t_rank_c Rank function for character @p c
//! \param t_cumulative_c Cumulative count for alphabet [0..sigma]
//! \param t_first Start position of queried range [t_first; t_last)
//! \param t_last End position of queried range [t_first; t_last)
//! \param t_c Queried character c
//! \param t_create_range Function to create new range corresponding to @p c in range [t_first; t_last)
//! \param t_empty_range Default empty range
//! \return New range corresponding to character @p c in range [t_first; t_last) or empty range if character does not appear in range
template<typename TRankC, typename TCumulativeC, typename TValue, typename TChar, typename TCreateNewRange, typename TNewRange>
auto lf(const TRankC &t_rank_c,
        const TCumulativeC &t_cumulative_c,
        const TValue &t_first,
        const TValue &t_last,
        const TChar &t_c,
        const TCreateNewRange &t_create_range,
        const TNewRange &t_empty_range) {
  // Number of c before the interval
  auto c_before_sp = t_rank_c(t_c, t_first);

  // Number of c before the interval + number of c inside the interval range
  auto c_until_ep = t_rank_c(t_c, t_last);

  // If there are no c in the interval, return empty range
  if (c_before_sp == c_until_ep)
    return t_empty_range;

  // Number of characters smaller than c
  auto smaller_c = t_cumulative_c(t_c);

  return t_create_range(c_before_sp, c_until_ep, smaller_c);
}

class EmptyClass {};

//! LF functor
//! \tparam TRankC Rank function for a given symbol
//! \tparam TCumulativeC Function to cumulative count for alphabet
//! \tparam TCreateRange Function to create range
//! \tparam TRange Range type
template<typename TRankC, typename TCumulativeC, typename TCreateRange, typename TRange>
class LF {
 public:
  LF(const TRankC &t_rank_c,
     const TCumulativeC &t_cumulative_c,
     const TCreateRange &t_create_range,
     const TRange &t_empty_range)
      : rank_c_{t_rank_c}, cumulative_c_{t_cumulative_c}, create_range_{t_create_range}, empty_range_{t_empty_range} {
  }

  //! LF function
  //! \tparam TTRange Range type {sp; ep}
  //! \tparam TChar Character type for compact alphabet
  //! \param t_range Range [sp; ep)
  //! \param t_c Character for compact alphabet in [0..sigma]
  //! \return [new_sp; new_ep)
  template<typename TTRange, typename TChar>
  auto operator()(const TTRange &t_range, const TChar &t_c) const {
    const auto&[first, last] = t_range;
    return (*this)(first, last, t_c);
  }

  //! LF function
  //! \tparam TValue Range position type {t_first; t_last}
  //! \tparam TChar Character type for compact alphabet
  //! \param t_first Start position of queried range [t_first; t_last)
  //! \param t_last End position of queried range [t_first; t_last)
  //! \param t_c Character for compact alphabet in [0..sigma]
  //! \return [new_sp; new_ep)
  template<typename TValue, typename TChar>
  auto operator()(const TValue &t_first, const TValue &t_last, const TChar &t_c) const {
    return lf(rank_c_, cumulative_c_, t_first, t_last, t_c, create_range_, empty_range_);
  }

 protected:
  TRankC rank_c_; // Rank function for partial psis
  TCumulativeC cumulative_c_; // Cumulative count for alphabet [0..sigma]

  TCreateRange create_range_; // Function to create new range
  TRange empty_range_; // Default new empty range
};

//! LF functor
//! \tparam TRankC Rank function for a given symbol
//! \tparam TCumulativeC Function to cumulative count for alphabet
template<typename TRankC, typename TCumulativeC>
class LF<TRankC, TCumulativeC, EmptyClass, EmptyClass> {
 public:
  LF(const TRankC &t_rank_c, const TCumulativeC &t_cumulative_c) : rank_c_{t_rank_c}, cumulative_c_{t_cumulative_c} {}

  //! LF function
  //! \tparam TTRange Range type {sp; ep}
  //! \tparam TChar Character type for compact alphabet
  //! \param t_range Range [sp; ep)
  //! \param t_c Character for compact alphabet in [0..sigma]
  //! \return [new_sp; new_ep)
  template<typename TTRange, typename TChar>
  auto operator()(const TTRange &t_range, const TChar &t_c) const {
    auto create_range = [](const auto &tt_c_before_sp, const auto &tt_c_until_ep, const auto &tt_smaller_c) {
      return TTRange{tt_smaller_c + tt_c_before_sp, tt_smaller_c + tt_c_until_ep};
    };
    TTRange empty_range{1, 0};

    const auto&[first, last] = t_range;
    return lf(rank_c_, cumulative_c_, first, last, t_c, create_range, empty_range);
  }

 protected:
  TRankC rank_c_; // Rank function for partial psis
  TCumulativeC cumulative_c_; // Cumulative count for alphabet [0..sigma]
};

//! LFOnPsi function using partial psi function
//! \tparam TRankPartialPsi Rank function for partial psis
//! \tparam TCumulativeC Function to cumulative count for alphabet
template<typename TRankPartialPsi, typename TCumulativeC>
class LFOnPsi : public LF<TRankPartialPsi, TCumulativeC, EmptyClass, EmptyClass> {
 public:
  LFOnPsi(const TRankPartialPsi &t_rank_partial_psi, const TCumulativeC &t_cumulative_c)
      : LF<TRankPartialPsi, TCumulativeC, EmptyClass, EmptyClass>{t_rank_partial_psi, t_cumulative_c} {
  }
};

}

#endif //SRI_LF_H_
