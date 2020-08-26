//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#ifndef RI_PREDECESSOR_H_
#define RI_PREDECESSOR_H_

#include <memory>
#include <cstddef>

#include <sdsl/rank_support.hpp>

namespace ri {

template<typename TBitvector, typename TRank = typename TBitvector::rank_1_type, typename TSelect = typename TBitvector::select_1_type>
class CircularPredecessor {
 public:
  CircularPredecessor(const std::shared_ptr<TBitvector> &t_bitvector,
                      const std::shared_ptr<TRank> &t_rank,
                      const std::shared_ptr<TSelect> &t_select)
      : bitvector_{t_bitvector}, rank_{t_rank}, select_{t_select} {
    n_ones_ = (*rank_)(t_bitvector->size());
  }

  /// Find the greatest number that is strictly minor than the value given (or the greatest number if i is lesser than
  /// all the numbers).
  /// \param i Number to search.
  /// \return <o, v>: the predecessor of i is the o-th number and its value is v.
  /// \note The search is circular. The predecessor is strictly lower (or higher if i is lesser all values).
  auto operator()(std::size_t i) const {
//    auto r = (*rank_)(i);
    auto order = ((*rank_)(i) + n_ones_ - 1) % n_ones_;
    auto value = (*select_)(order + 1);

    return std::make_pair(order, value);
  }

 private:
  std::shared_ptr<TBitvector> bitvector_;
  std::shared_ptr<TRank> rank_;
  std::shared_ptr<TSelect> select_;

  std::size_t n_ones_;
};

}

#endif //RI_PREDECESSOR_H_
