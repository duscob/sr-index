//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "sr-index/rle_string.hpp"

using BWT = std::string;
using Runs = std::vector<ri::StringRun>;

namespace ri {
void PrintTo(const StringRun &run, std::ostream *os) {
  *os << "{run: " << run.run << "; char: " << run.c << "; range: [" << run.range.first << "; " << run.range.second
      << "]}";
}
}

class RLEString_Tests : public testing::TestWithParam<std::tuple<BWT, ri::range_t, Runs>> {
};

TEST_P(RLEString_Tests, break_in_runs) {
  const auto &bwt = std::get<0>(GetParam());
  ri::rle_string<> bwt_rle(bwt);

  const auto &range = std::get<1>(GetParam());
  auto runs_in_range = bwt_rle.break_in_runs(range);

  const auto &e_runs_in_range = std::get<2>(GetParam());
  EXPECT_THAT(runs_in_range, testing::ElementsAreArray(e_runs_in_range));
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RLEString_Tests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{0, 11},
                        Runs{ri::StringRun{0, 4, {0, 1}}, ri::StringRun{1, 3, {2, 2}}, ri::StringRun{2, 4, {3, 3}},
                             ri::StringRun{3, 1, {4, 4}}, ri::StringRun{4, 2, {5, 8}}, ri::StringRun{5, 3, {9, 11}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{4, 11},
                        Runs{ri::StringRun{3, 1, {4, 4}}, ri::StringRun{4, 2, {5, 8}}, ri::StringRun{5, 3, {9, 11}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{1, 11},
                        Runs{ri::StringRun{0, 4, {1, 1}}, ri::StringRun{1, 3, {2, 2}}, ri::StringRun{2, 4, {3, 3}},
                             ri::StringRun{3, 1, {4, 4}}, ri::StringRun{4, 2, {5, 8}}, ri::StringRun{5, 3, {9, 11}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{1, 7},
                        Runs{ri::StringRun{0, 4, {1, 1}}, ri::StringRun{1, 3, {2, 2}}, ri::StringRun{2, 4, {3, 3}},
                             ri::StringRun{3, 1, {4, 4}}, ri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{5, 8},
                        Runs{ri::StringRun{4, 2, {5, 8}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{6, 8},
                        Runs{ri::StringRun{4, 2, {6, 8}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{5, 7},
                        Runs{ri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{6, 7},
                        Runs{ri::StringRun{4, 2, {6, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{7, 7},
                        Runs{ri::StringRun{4, 2, {7, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{3, 3},
                        Runs{ri::StringRun{2, 4, {3, 3}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{5, 10},
                        Runs{ri::StringRun{4, 2, {5, 8}}, ri::StringRun{5, 3, {9, 10}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        ri::range_t{7, 9},
                        Runs{ri::StringRun{4, 2, {7, 8}}, ri::StringRun{5, 3, {9, 9}}})
    )
);
