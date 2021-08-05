//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/9/20.
//

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sr-index/bwt.h"
#include "sr-index/rle_string.hpp"

using Sequence = std::string;
using SA = std::vector<std::size_t>;
using BWT = std::string;
using BWTRunHeads = std::vector<std::tuple<std::size_t, std::size_t, char, std::size_t>>;
using BWTRunTails = std::vector<std::tuple<std::size_t, std::size_t, char, std::size_t>>;

class ComputeBWT_Tests : public testing::TestWithParam<std::tuple<Sequence, SA, BWT, BWTRunHeads, BWTRunTails>> {
};

TEST_P(ComputeBWT_Tests, FromSAAndText) {
  const auto &seq = std::get<0>(GetParam());
  const auto &sa = std::get<1>(GetParam());

  BWT bwt;
  auto report_bwt = [&bwt](auto idx, auto symbol) { bwt.push_back(symbol); };
  BWTRunHeads bwt_heads;
  auto report_bwt_head = [&bwt_heads](const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
    bwt_heads.emplace_back(run, idx, symbol, pos);
  };
  BWTRunTails bwt_tails;
  auto report_bwt_tail = [&bwt_tails](const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
    bwt_tails.emplace_back(run, idx, symbol, pos);
  };

  auto n_runs = sri::computeBWT(seq.size(),
                                [&sa](auto idx) { return sa[idx]; },
                                [&seq](auto idx) { return seq[idx]; },
                                report_bwt,
                                report_bwt_head,
                                report_bwt_tail);

  const auto &e_bwt = std::get<2>(GetParam());
  const auto &e_bwt_heads = std::get<3>(GetParam());
  const auto &e_bwt_tails = std::get<4>(GetParam());

  EXPECT_EQ(n_runs, e_bwt_heads.size());
  EXPECT_EQ(bwt, e_bwt);
  EXPECT_THAT(bwt_heads, testing::ElementsAreArray(e_bwt_heads));
  EXPECT_THAT(bwt_tails, testing::ElementsAreArray(e_bwt_tails));
}

INSTANTIATE_TEST_SUITE_P(
    BWT,
    ComputeBWT_Tests,
    testing::Values(
        std::make_tuple("abracadabra$",
                        SA{11, 10, 7, 0, 3, 5, 8, 1, 4, 6, 9, 2},
                        "ard$rcaaaabb",
                        BWTRunHeads{{0, 0, 'a', 10}, {1, 1, 'r', 9}, {2, 2, 'd', 6}, {3, 3, '$', 11}, {4, 4, 'r', 2},
                                    {5, 5, 'c', 4}, {6, 6, 'a', 7}, {7, 10, 'b', 8}},
                        BWTRunTails{{0, 0, 'a', 10}, {1, 1, 'r', 9}, {2, 2, 'd', 6}, {3, 3, '$', 11}, {4, 4, 'r', 2},
                                    {5, 5, 'c', 4}, {6, 9, 'a', 5}, {7, 11, 'b', 1}}
        ),
        std::make_tuple("abcabcabcab$",
                        SA{11, 9, 6, 3, 0, 10, 7, 4, 1, 8, 5, 2},
                        "bccc$aaaabbb",
                        BWTRunHeads{{0, 0, 'b', 10}, {1, 1, 'c', 8}, {2, 4, '$', 11}, {3, 5, 'a', 9}, {4, 9, 'b', 7}},
                        BWTRunTails{{0, 0, 'b', 10}, {1, 3, 'c', 2}, {2, 4, '$', 11}, {3, 8, 'a', 0}, {4, 11, 'b', 1}}
        )
    )
);


using F = std::vector<std::size_t>;
using Range = std::pair<std::size_t, std::size_t>;
using Char = unsigned char;

class LFOnBWT_Tests : public testing::TestWithParam<std::tuple<BWT, F, Range, Char, Range>> {
};

TEST_P(LFOnBWT_Tests, compute) {
  const auto &bwt = std::get<0>(GetParam());

  sri::rle_string<> bwt_rle(bwt);
  auto get_bwt_rank = [&bwt_rle](auto pos, auto c) {
    return bwt_rle.rank(pos, c);
  };

  const auto &f = std::get<1>(GetParam());
  auto get_f = [&f](auto i) {
    return f[i];
  };

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());
//  Char max_c = *std::max_element(bwt.begin(), bwt.end()) + 1;
  Char max_c = bwt.size() - 1;
  auto new_range = sri::computeLF(get_bwt_rank, get_f, range, c, bwt.size(), max_c);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

INSTANTIATE_TEST_SUITE_P(
    BWT,
    LFOnBWT_Tests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{0, 11}, 0, Range{1, 0}), // 0 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{0, 11}, 1, Range{0, 0}), // 1 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{0, 11}, 2, Range{1, 4}), // 2 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{0, 11}, 3, Range{5, 8}), // 3 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{0, 11}, 4, Range{9, 11}), // 4 before empty

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{1, 4}, 1, Range{0, 0}), // 1 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{1, 4}, 2, Range{1, 0}), // 2 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{1, 4}, 3, Range{5, 5}), // 3 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{1, 4}, 4, Range{10, 11}), // 4 before 2

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{10, 11}, 2, Range{1, 0}), // 2 before 42
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, F{0, 0, 1, 5, 9, 12}, Range{10, 11}, 3, Range{7, 8}) // 3 before 42
    )
);
