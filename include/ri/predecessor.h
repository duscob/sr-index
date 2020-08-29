//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#ifndef RI_PREDECESSOR_H_
#define RI_PREDECESSOR_H_

#include <memory>
#include <cstddef>

#include <sdsl/rank_support.hpp>

namespace ri {

template<typename TRank, typename TSelect>
class CircularPredecessor {
 public:
  CircularPredecessor(const TRank &t_rank, const TSelect &t_select, std::size_t t_bitvector_size)
      : rank_{t_rank}, select_{t_select} {
    n_ones_ = rank_(t_bitvector_size);
  }

  /// Find the greatest number that is strictly minor than the value given (or the greatest number if i is lesser than
  /// all the numbers).
  /// \param i Number to search.
  /// \return <o, v>: the predecessor of i is the o-th number and its value is v.
  /// \note The search is circular. The predecessor is strictly lower (or higher if i is lesser all values).
  auto operator()(std::size_t i) const {
    auto order = (rank_(i) + n_ones_ - 1) % n_ones_;
    auto value = select_(order + 1);

    return std::make_pair(order, value);
  }

 private:
  TRank rank_;
  TSelect select_;

  std::size_t n_ones_;
};

template<typename TRank, typename TSelect>
auto buildCircularPredecessor(const TRank & t_rank, const TSelect &t_select, std::size_t t_bitvector_size) {
  return CircularPredecessor<TRank, TSelect>(t_rank, t_select, t_bitvector_size);
}

}

#endif //RI_PREDECESSOR_H_
