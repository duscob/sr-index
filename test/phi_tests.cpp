//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

#include "sr-index/phi.h"
#include "sr-index/sequence_ops.h"
#include "sr-index/tools.h"
#include "sr-index/rle_string.hpp"
#include "sr-index/bwt.h"

using BitVector = sdsl::bit_vector;
using IntVector = sdsl::int_vector<>;
using Param = std::size_t;
using Result = std::pair<std::size_t, bool>;
using ParamResult = std::pair<Param, Result>;

class Phi_Tests : public testing::TestWithParam<std::tuple<BitVector, IntVector, BitVector, IntVector, ParamResult>> {
};

class PhiBackward_Tests : public Phi_Tests {
};

TEST_P(PhiBackward_Tests, compute) {
  const auto &bv = std::get<0>(GetParam());
  auto rank = sdsl::bit_vector::rank_1_type(&bv);
  auto select = sdsl::bit_vector::select_1_type(&bv);
  auto predecessor = sri::buildCircularPredecessor(std::ref(rank), std::ref(select), bv.size());

  const auto &pred_to_run = std::get<1>(GetParam());
  const auto &is_trustworthy = std::get<2>(GetParam());
  auto get_pred_to_run = sri::buildRandomAccessForTwoContainers(std::ref(pred_to_run), std::ref(is_trustworthy));

  const auto &samples = std::get<3>(GetParam());
  auto get_sample = sri::buildRandomAccessForContainer(std::ref(samples));
  sri::SampleValidatorDefault sampled_tail_validator_default;

  auto phi = sri::buildPhiBackward(predecessor, get_pred_to_run, get_sample, sampled_tail_validator_default, bv.size());

//  auto phi = BuildPhi(GetParam());

  const auto &prev_value = std::get<4>(GetParam()).first;
  auto value = phi(prev_value);

  const auto &e_value = std::get<4>(GetParam()).second;
  EXPECT_EQ(value, e_value);
}

INSTANTIATE_TEST_SUITE_P(
    WithoutSampling,
    PhiBackward_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1}), // BWT heads
        testing::Values(IntVector{1, 3, 0, 4, 5, 2}), // Predecessor to sample
        testing::Values(BitVector{1, 1, 1, 1, 1, 1}), // Trustworthy sample?
        testing::Values(IntVector{5, 7, 2, 11, 0, 1}), // Sampled tails
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
    WithSampling,
    PhiBackward_Tests,
    testing::Combine(
        testing::Values(BitVector{0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1}), // BWT heads
        testing::Values(IntVector{0, 2, 1}), // Predecessor to sample
        testing::Values(BitVector{1, 0, 1}), // Trustworthy sample?
        testing::Values(IntVector{7, 2, 11}), // Sampled tails
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
//    PhiInvWithoutSampling,
//    PhiBackward_Tests,
//    testing::Combine(
//        testing::Values(BitVector{1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1}), // BWT tails
//        testing::Values(IntVector{5, 0, 3, 1, 2, 4}), // Predecessor to sample
//        testing::Values(BitVector{1, 1, 1, 1, 1, 1}), // Trustworthy sample?
//        testing::Values(IntVector{10, 7, 2, 11, 6, 9}), // Sampled heads
//        testing::Values(//ParamResult{11, {6, true}},
//            ParamResult{6, {8, true}},
//            ParamResult{8, {3, true}},
//            ParamResult{3, {0, true}},
//            ParamResult{0, {7, true}},
//            ParamResult{7, {9, true}},
//            ParamResult{9, {4, true}},
//            ParamResult{4, {1, true}},
//            ParamResult{1, {10, true}},
//            ParamResult{10, {5, true}},
//            ParamResult{5, {2, true}}
//        )
//    )
//);

class PhiForward_Tests : public Phi_Tests {
};

TEST_P(PhiForward_Tests, compute) {
  const auto &bv = std::get<0>(GetParam());
  auto rank = sdsl::bit_vector::rank_1_type(&bv);
  auto select = sdsl::bit_vector::select_1_type(&bv);
  auto successor = sri::CircularSoftSuccessor(std::ref(rank), std::ref(select), bv.size());

  const auto &sample_index = std::get<1>(GetParam());
  const auto &is_trustworthy = std::get<2>(GetParam());
  auto get_sample_index = sri::buildRandomAccessForTwoContainers(std::ref(sample_index), std::ref(is_trustworthy));

  const auto &samples = std::get<3>(GetParam());
  auto get_sample = sri::buildRandomAccessForContainer(std::ref(samples));
  sri::SampleValidatorDefault sampled_tail_validator_default;

  auto phi = sri::buildPhiForward(successor, get_sample_index, get_sample, sampled_tail_validator_default, bv.size());

  const auto &prev_value = std::get<4>(GetParam()).first;
  auto value = phi(prev_value);

  const auto &e_value = std::get<4>(GetParam()).second;
  EXPECT_EQ(value, e_value);
}

INSTANTIATE_TEST_SUITE_P(
    WithoutSampling,
    PhiForward_Tests,
    testing::Combine(
        testing::Values(BitVector{1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1}), // BWT tails
        testing::Values(IntVector{1, 0, 3, 2, 5, 4}), // Successor to sample
        testing::Values(BitVector{1, 1, 1, 1, 1, 1}), // Trustworthy sample?
        testing::Values(IntVector{10, 7, 2, 11, 6, 9}), // Sampled heads
        testing::Values(ParamResult{11, {6, true}},
                        ParamResult{6, {8, true}},
                        ParamResult{8, {3, true}},
                        ParamResult{3, {0, true}},
                        ParamResult{0, {7, true}},
                        ParamResult{7, {9, true}},
                        ParamResult{9, {4, true}},
                        ParamResult{4, {1, true}},
                        ParamResult{1, {10, true}},
                        ParamResult{10, {5, true}},
                        ParamResult{5, {2, true}}
        )
    )
);

INSTANTIATE_TEST_SUITE_P(
    WithSampling,
    PhiForward_Tests,
    testing::Combine(
        testing::Values(BitVector{1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1}), // BWT tails
        testing::Values(IntVector{0, 2, 1, 3}), // Successor to sample
        testing::Values(BitVector{1, 0, 1, 0}), // Trustworthy sample?
        testing::Values(IntVector{7, 2, 11, 6}), // Sub-sampled heads  (sub-sampled indexes = {3, 1, 4, 2})
        testing::Values(ParamResult{11, {6, false}},
                        ParamResult{6, {1, false}},
                        ParamResult{8, {3, false}},
                        ParamResult{3, {0, true}},
                        ParamResult{0, {7, true}},
                        ParamResult{7, {2, false}},
                        ParamResult{9, {4, false}},
                        ParamResult{4, {1, true}},
                        ParamResult{1, {10, false}},
                        ParamResult{10, {5, false}},
                        ParamResult{5, {2, true}}
        )
    )
);

using BWT = std::string;
using F = sdsl::int_vector<>;
using Psi = IntVector;
using BWTHeadsPos = BitVector;
using BWTTailsPos = BitVector;
using Links = IntVector;

class ComputeLinkForPhiForward_Tests
    : public testing::TestWithParam<std::tuple<BWT, F, Psi, std::size_t, BWTHeadsPos, BWTTailsPos, Links>> {
};

TEST_P(ComputeLinkForPhiForward_Tests, MarkToSample_bv) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  const auto &f = std::get<1>(GetParam());
  auto get_f = sri::buildRandomAccessForContainer(std::cref(f));
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  auto psi = sri::RandomAccessForCRefContainer(std::cref(std::get<2>(GetParam())));

  auto r = std::get<3>(GetParam());

  const auto &bwt_heads_pos = std::get<4>(GetParam());
  auto heads_rank = BWTHeadsPos::rank_1_type(&bwt_heads_pos);

  const auto &bwt_tails_pos = std::get<5>(GetParam());
  auto tails_select = BWTTailsPos::select_1_type(&bwt_tails_pos);

  IntVector links(r, 0);
  auto report_link = [&links](auto tt_i, auto tt_l) { links[tt_i] = tt_l; };

  sri::computeMarkToSampleLinkForPhiForward(bwt.size(), r, tails_select, heads_rank, lf, psi, report_link);

  const auto &e_links = std::get<6>(GetParam());
  EXPECT_THAT(links, testing::ElementsAreArray(e_links));
}

TEST_P(ComputeLinkForPhiForward_Tests, MarkToSample_v) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  const auto &f = std::get<1>(GetParam());
  auto get_f = sri::buildRandomAccessForContainer(std::cref(f));
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  auto psi = sri::RandomAccessForCRefContainer(std::cref(std::get<2>(GetParam())));

  auto r = std::get<3>(GetParam());

  auto get_mark_pos = [](const auto &tt_bv) {
    std::vector<std::size_t> pos;
    for (int i = 0; i < tt_bv.size(); ++i) {
      if (tt_bv[i])
        pos.push_back(i);
    }
    return pos;
  };

  auto bwt_head_pos = get_mark_pos(std::get<4>(GetParam()));
  auto rank_sample = [&bwt_head_pos](const auto &tt_k) {
    return std::lower_bound(bwt_head_pos.begin(), bwt_head_pos.end(), tt_k) - bwt_head_pos.begin();
  };

  auto bwt_tail_pos = get_mark_pos(std::get<5>(GetParam()));

  IntVector links(r, 0);
  for (std::size_t i = 0; i < r; ++i) {
    links[i] = sri::computeMarkToSampleLinkForPhiForward(bwt_tail_pos[i], bwt.size(), lf, psi, rank_sample);
  }

  const auto &e_links = std::get<6>(GetParam());
  EXPECT_THAT(links, testing::ElementsAreArray(e_links));
}

TEST_P(ComputeLinkForPhiForward_Tests, SampleToMark_v) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  const auto &f = std::get<1>(GetParam());
  auto get_f = sri::buildRandomAccessForContainer(std::cref(f));
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  auto psi = sri::RandomAccessForCRefContainer(std::cref(std::get<2>(GetParam())));

  auto r = std::get<3>(GetParam());

  auto get_mark_pos = [](const auto &tt_bv) {
    std::vector<std::size_t> pos;
    for (int i = 0; i < tt_bv.size(); ++i) {
      if (tt_bv[i])
        pos.push_back(i);
    }
    return pos;
  };

  auto bwt_head_pos = get_mark_pos(std::get<4>(GetParam()));

  auto e_bwt_tail_pos = get_mark_pos(std::get<5>(GetParam()));
  auto rank_mark = [&e_bwt_tail_pos](const auto &tt_k) {
    return std::lower_bound(e_bwt_tail_pos.begin(), e_bwt_tail_pos.end(), tt_k) - e_bwt_tail_pos.begin();
  };

  const auto &links = std::get<6>(GetParam());

  IntVector bwt_tail_pos(r, 0);
  for (std::size_t i = 0; i < r; ++i) {
    bwt_tail_pos[i] = sri::computeSampleToMarkLinkForPhiForward(bwt_head_pos[links[i]], bwt.size(), lf, psi, rank_mark);
  }

  EXPECT_EQ(bwt_tail_pos.size(), e_bwt_tail_pos.size());
  for (int i = 0; i < bwt_tail_pos.size(); ++i) {
    EXPECT_EQ(bwt_tail_pos[i], rank_mark(e_bwt_tail_pos[i])) << "error at position " << i;
  }
}

INSTANTIATE_TEST_SUITE_P(
    Links,
    ComputeLinkForPhiForward_Tests,
    testing::Values(
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
            F{0, 0, 1, 5, 9, 12},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
            6,
            BWTHeadsPos{1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0},
            BWTTailsPos{0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1},
            Links{2, 5, 3, 4, 1, 0}
        )
    )
);

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
  // PhiBackward
  const auto &bv = std::get<0>(GetParam());
  auto rank = sdsl::bit_vector::rank_1_type(&bv);
  auto select = sdsl::bit_vector::select_1_type(&bv);
  auto predecessor = sri::buildCircularPredecessor(std::ref(rank), std::ref(select), bv.size());

  const auto &pred_to_run = std::get<1>(GetParam());
  const auto &is_trustworthy = std::get<2>(GetParam());
  auto get_pred_to_run = sri::buildRandomAccessForTwoContainers(std::ref(pred_to_run), std::ref(is_trustworthy));

  const auto &samples = std::get<3>(GetParam());
  auto get_sample = sri::buildRandomAccessForContainer(std::ref(samples));
  sri::SampleValidatorDefault sampled_tail_validator_default;

  auto phi = sri::buildPhiBackward(predecessor, get_pred_to_run, get_sample, sampled_tail_validator_default, bv.size());

  // Split in runs
  const auto &bwt = std::get<4>(GetParam());
  sri::rle_string<> bwt_rle(bwt);
  auto split_in_runs = sri::buildSplitInRuns(std::ref(bwt_rle));

  // Backward navigation
  auto get_rank_of_char = sri::buildRankOfChar(std::ref(bwt_rle));
  const auto &f = std::get<5>(GetParam());
  auto get_f = sri::buildRandomAccessForContainer(std::ref(f));
  auto max_c = f.size() - 1;
  auto lf = sri::buildLF(get_rank_of_char, std::ref(get_f), bwt.size(), max_c);

  // Get sample
  const auto &is_run_sampled = std::get<6>(GetParam());
  BitVector::rank_1_type rank_run_tail(&is_run_sampled);
  auto get_sample_for_bwt_run = sri::buildGetSampleForBWTRun(std::ref(rank_run_tail), get_sample);

  auto run_of_sa_pos = sri::buildRunOfSAPosition(std::ref(bwt_rle));
  auto get_is_run_sampled = sri::buildRandomAccessForContainer(std::ref(is_run_sampled));
  auto get_sample_for_sa_pos =
      sri::buildGetSampleForSAPosition(std::ref(run_of_sa_pos), get_is_run_sampled, get_sample_for_bwt_run);

  const auto &sampling_size = std::get<7>(GetParam());

  auto phi_for_range = sri::buildPhiForRange(phi, split_in_runs, lf, get_sample_for_sa_pos, sampling_size, bwt.size());

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

        testing::Values(5),

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