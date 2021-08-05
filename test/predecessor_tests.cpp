//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#include <algorithm>

#include <gtest/gtest.h>

#include <sdsl/bit_vectors.hpp>

#include "sr-index/predecessor.h"

using BitVector = std::vector<bool>;
using Value = std::size_t;
using Predecessor = std::pair<std::size_t, std::size_t>;

class CircularPredecessor_Tests : public testing::TestWithParam<std::tuple<BitVector, Value, Predecessor>> {
};

TEST_P(CircularPredecessor_Tests, compute) {
  const auto &raw_bv = std::get<0>(GetParam());

  auto bv = sdsl::bit_vector(raw_bv.size());
  copy(raw_bv.begin(), raw_bv.end(), bv.begin());

  auto rank = sdsl::bit_vector::rank_1_type(&bv);
  auto select = sdsl::bit_vector::select_1_type(&bv);

  auto predecessor = ri::buildCircularPredecessor(std::ref(rank), std::ref(select), bv.size());

  const auto &value = std::get<1>(GetParam());
  auto p = predecessor(value);

  const auto &e_p = std::get<2>(GetParam());
  EXPECT_EQ(e_p, p);
}

INSTANTIATE_TEST_SUITE_P(
    Predecessor,
    CircularPredecessor_Tests,
    testing::Values(
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 0, Predecessor{5, 11}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 1, Predecessor{5, 11}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 2, Predecessor{5, 11}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 3, Predecessor{0, 2}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 4, Predecessor{0, 2}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 5, Predecessor{0, 2}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 6, Predecessor{0, 2}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 7, Predecessor{1, 6}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 8, Predecessor{2, 7}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 9, Predecessor{2, 7}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 10, Predecessor{3, 9}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 11, Predecessor{4, 10}),
        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}, 12, Predecessor{5, 11})
    )
);


