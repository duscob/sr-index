//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#ifndef SRI_SEQUENCE_OPS_H_
#define SRI_SEQUENCE_OPS_H_

#include <memory>
#include <cstddef>

namespace sri {

//! Find the greatest number strictly smaller than the value i or the greatest number if i is smaller than all number.
/**
 * @tparam TRank Rank function for the bit-vector with the marked values
 * @tparam TSelect Select function for the bit-vector with the marked values
 */
template<typename TRank, typename TSelect>
class CircularPredecessor {
 public:
  CircularPredecessor(const TRank &t_rank, const TSelect &t_select, std::size_t t_bitvector_size)
      : rank_{t_rank}, select_{t_select} {
    n_ones_ = rank_(t_bitvector_size);
  }

  //! Predecessor function
  /**
   * @param i Number to search.
   * @return <o, v>: the predecessor of i is the o-th number and its value is v.
   * @note The search is circular. The predecessor is strictly lower (or higher if i is less than all values).
   */
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
auto buildCircularPredecessor(const TRank &t_rank, const TSelect &t_select, std::size_t t_bitvector_size) {
  return CircularPredecessor<TRank, TSelect>(t_rank, t_select, t_bitvector_size);
}

//! Find the smallest number greater or equal than the value i or the smallest number if i is greater than all numbers.
/**
 * @tparam TRank Rank function for the bit-vector with the marked values
 * @tparam TSelect Select function for the bit-vector with the marked values
 */
template<typename TRank, typename TSelect>
class CircularSoftSuccessor {
 public:
  CircularSoftSuccessor(const TRank &t_rank, const TSelect &t_select, std::size_t t_bitvector_size)
      : rank_{t_rank}, select_{t_select} {
    n_ones_ = rank_(t_bitvector_size);
  }

  //! Soft successor function
  /**
   * @param i Number to search.
   * @return <o, v>: the successor of i is the o-th number and its value is v.
   * @note The search is circular. The soft successor is higher or equal (or lower if i is greater all values).
   */
  auto operator()(std::size_t i) const {
    auto order = rank_(i) % n_ones_;
    auto value = select_(order + 1);

    return std::make_pair(order, value);
  }

 private:
  TRank rank_;
  TSelect select_;

  std::size_t n_ones_;
};
}

#endif //SRI_SEQUENCE_OPS_H_
