//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/bit_vectors.hpp>

#include <ri/phi.h>
#include <ri/predecessor.h>
#include <ri/tools.h>
#include <ri/rle_string.hpp>
#include <ri/bwt.h>

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

//  auto phi = BuildPhi(GetParam());

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

using BWT = std::string;
using F = sdsl::int_vector<>;
using SamplingSize = std::size_t;
using Range = std::pair<std::size_t, std::size_t>;
using RangeXPrevValue = std::pair<Range, std::size_t>;
using SA = std::vector<std::size_t>;
using RangeXPrevValueXSA = std::pair<RangeXPrevValue, SA>;

class PhiForRange_Tests
    : public testing::TestWithParam<std::tuple<BitVector,
                                               IntVector,
                                               BitVector,
                                               IntVector,
                                               BWT,
                                               F,
                                               BitVector,
                                               SamplingSize,
                                               RangeXPrevValueXSA>> {
};

TEST_P(PhiForRange_Tests, compute) {
  // Phi
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

  // Split in runs
  const auto &bwt = std::get<4>(GetParam());
  ri::rle_string<> bwt_rle(bwt);
  auto split_in_runs = ri::buildSplitInRuns(std::ref(bwt_rle));

  // Backward navigation
  auto get_rank_of_char = ri::buildRankOfChar(std::ref(bwt_rle));
  const auto &f = std::get<5>(GetParam());
  auto get_f = ri::buildRandomAccessForContainer(std::ref(f));
  auto max_c = f.size() - 1;
  auto lf = ri::buildLF(get_rank_of_char, std::ref(get_f), bwt.size(), max_c);

  // Get sample
  const auto &is_run_sampled = std::get<6>(GetParam());
  BitVector::rank_1_type rank_run_tail(&is_run_sampled);
  auto get_sample_for_bwt_run = ri::buildGetSampleForBWTRun(std::ref(rank_run_tail), get_sample);

  auto run_of_sa_pos = ri::buildRunOfSAPosition(std::ref(bwt_rle));
  auto get_is_run_sampled = ri::buildRandomAccessForContainer(std::ref(is_run_sampled));
  auto get_sample_for_sa_pos =
      ri::buildGetSampleForSAPosition(std::ref(run_of_sa_pos), get_is_run_sampled, get_sample_for_bwt_run);

  const auto &sampling_size = std::get<7>(GetParam());

  auto phi_for_range = ri::buildPhiForRange(phi, split_in_runs, lf, get_sample_for_sa_pos, sampling_size, bwt.size());

  const auto &input = std::get<8>(GetParam()).first;
  std::vector<std::size_t> sa;
  auto report = [&sa](const auto &value) { sa.emplace_back(value); };
  phi_for_range(input.first, input.second, report);

  const auto &e_sa = std::get<8>(GetParam()).second;
  EXPECT_THAT(sa, testing::ElementsAreArray(e_sa));
}

INSTANTIATE_TEST_SUITE_P(
    PhiForRangeWithoutSampling,
    PhiForRange_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}),
        testing::Values(IntVector{1, 3, 0, 4, 5, 2}),
        testing::Values(BitVector{1, 1, 1, 1, 1, 1}),
        testing::Values(IntVector{5, 7, 2, 11, 0, 1}),

        testing::Values(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}),
        testing::Values(F{0, 0, 1, 5, 9, 12}),

        testing::Values(BitVector{1, 1, 1, 1, 1, 1}),

        testing::Values(0),

        testing::Values(
            RangeXPrevValueXSA{RangeXPrevValue{Range{1, 10}, 2}, SA{5, 10, 1, 4, 9, 7, 0, 3, 8, 6}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{5, 9}, 5}, SA{10, 1, 4, 9, 7}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{2, 3}, 0}, SA{3, 8}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{1, 1}, 8}, SA{6}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{4, 4}, 7}, SA{0}}
        )
    )
);

INSTANTIATE_TEST_SUITE_P(
    PhiForRangeWithSampling,
    PhiForRange_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}),
        testing::Values(IntVector{0, 2, 1}),
        testing::Values(BitVector{1, 0, 1}),
        testing::Values(IntVector{7, 2, 11}),

        testing::Values(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}),
        testing::Values(F{0, 0, 1, 5, 9, 12}),

        testing::Values(BitVector{0, 1, 1, 1, 0, 0}),

        testing::Values(4),

        testing::Values(
            RangeXPrevValueXSA{RangeXPrevValue{Range{1, 10}, 2}, SA{5, 10, 1, 4, 9, 7, 0, 3, 8, 6}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{5, 9}, 5}, SA{10, 1, 4, 9, 7}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{2, 3}, 0}, SA{3, 8}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{1, 1}, 8}, SA{6}},
            RangeXPrevValueXSA{RangeXPrevValue{Range{4, 4}, 7}, SA{0}}
        )
    )
);
//SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}
//IntVector{2, 4, 1, 5, 0, 3}