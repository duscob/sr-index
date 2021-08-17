//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/17/21.
//

#include <string>

#include <gtest/gtest.h>

#include "sr-index/toehold.h"
#include "sr-index/rle_string.hpp"

using BWT = std::string;
using Range = std::pair<std::size_t, std::size_t>;
using Char = char;
using Step = std::size_t;
using Data = sri::DataBackwardSearchStep<Char>;

auto tie(const Data &t_data) {
  return std::tie(t_data.c, t_data.step, t_data.range_start, t_data.range_end);
}

class GetLastSpecialBackwardSearchStep_Test
    : public testing::TestWithParam<std::tuple<BWT, std::tuple<Range, Range, Char, Step, Data>, Data>> {
};

TEST_P(GetLastSpecialBackwardSearchStep_Test, execute) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &[range, next_range, c, step, prev_data] = std::get<1>(GetParam());

  auto get_data = sri::buildGetLastSpecialBackwardSearchStepForPhiBackward(std::cref(bwt));
  auto data = get_data(range, next_range, c, step, prev_data);

  const auto &e_data = std::get<2>(GetParam());
  EXPECT_EQ(tie(data), tie(e_data));
}

INSTANTIATE_TEST_SUITE_P(
    Toehold,
    GetLastSpecialBackwardSearchStep_Test,
    testing::Values(
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{0, 0}, 1, 0, Data{3, 0, 0, 11}), // Input
            Data{1, 0, 0, 11} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{1, 5}, 2, 0, Data{3, 0, 0, 11}), // Input
            Data{2, 0, 0, 11} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{9, 11}, 4, 0, Data{3, 0, 0, 11}), // Input
            Data{4, 0, 0, 11} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{9, 11}, 4, 2, Data{3, 2, 0, 11}), // Input
            Data{4, 2, 0, 11} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{9, 11}, Range{6, 8}, 3, 1, Data{4, 2, 0, 11}), // Input
            Data{4, 2, 0, 11} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{9, 11}, Range{1, 0}, 1, 1, Data{4, 2, 0, 11}), // Input
            Data{4, 2, 0, 11} // Next data
        )
    )
);