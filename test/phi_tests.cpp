//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#include <gtest/gtest.h>

#include <sdsl/bit_vectors.hpp>

#include <ri/phi.h>
#include <ri/predecessor.h>
#include <ri/tools.h>

using BitVector = sdsl::bit_vector;
using IntVector = sdsl::int_vector<>;
using Param = std::size_t;
using Result = std::pair<std::size_t, bool>;
using ParamResult = std::pair<Param, Result>;

class Phi_Tests : public testing::TestWithParam<std::tuple<BitVector, IntVector, BitVector, IntVector, ParamResult>> {
};

TEST_P(Phi_Tests, compute) {
  const auto &bv = std::get<0>(GetParam());
  auto rank = sdsl::bit_vector::rank_1_type(&bv);
  auto select = sdsl::bit_vector::select_1_type(&bv);
  auto predecessor = ri::buildCircularPredecessor(std::ref(rank), std::ref(select), bv.size());

  const auto &pred_to_run = std::get<1>(GetParam());
  const auto &is_trustworthy = std::get<2>(GetParam());
  auto get_pred_to_run = ri::buildRandomAccessForTwoContainers(std::ref(pred_to_run), std::ref(is_trustworthy));

  const auto &samples = std::get<3>(GetParam());
  auto get_sample = ri::buildRandomAccessForContainer(std::ref(samples));

  auto phi = ri::buildPhi(predecessor, get_pred_to_run, get_sample, bv.size());

  const auto &prev_value = std::get<4>(GetParam()).first;
  auto value = phi(prev_value);

  const auto &e_value = std::get<4>(GetParam()).second;
  EXPECT_EQ(value, e_value);
}

INSTANTIATE_TEST_SUITE_P(
    PhiWithoutSampling,
    Phi_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}),
        testing::Values(IntVector{1, 3, 0, 4, 5, 2}),
        testing::Values(BitVector{1, 1, 1, 1, 1, 1}),
        testing::Values(IntVector{5, 7, 2, 11, 0, 1}),
        testing::Values(ParamResult{2, {5, true}},
                        ParamResult{5, {10, true}},
                        ParamResult{10, {1, true}},
                        ParamResult{1, {4, true}},
                        ParamResult{4, {9, true}},
                        ParamResult{9, {7, true}},
                        ParamResult{7, {0, true}},
                        ParamResult{0, {3, true}},
                        ParamResult{3, {8, true}},
                        ParamResult{8, {6, true}},
                        ParamResult{6, {11, true}})
    )
);

INSTANTIATE_TEST_SUITE_P(
    PhiWithSampling,
    Phi_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}),
        testing::Values(IntVector{0, 2, 1}),
        testing::Values(BitVector{1, 0, 1}),
        testing::Values(IntVector{7, 2, 11}),
        testing::Values(ParamResult{2, {5, true}},
                        ParamResult{5, {10, true}},
                        ParamResult{10, {3, false}},
                        ParamResult{1, {4, true}},
                        ParamResult{4, {9, true}},
                        ParamResult{9, {2, false}},
                        ParamResult{7, {0, false}},
                        ParamResult{0, {3, true}},
                        ParamResult{3, {8, true}},
                        ParamResult{8, {1, false}},
                        ParamResult{6, {11, true}})
    )
);
//INSTANTIATE_TEST_SUITE_P(
//    Phi,
//    Phi_Tests,
//    testing::Values(
//        std::make_tuple(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1},
//                        IntVector{1, 3, 0, 4, 5, 2},
//                        BitVector{1, 1, 1, 1, 1, 1},
//                        IntVector{5, 7, 2, 11, 0, 1},
//                        ParamResult{6, {11, true}})
//    )
//);

//IntVector{2, 4, 1, 5, 0, 3}