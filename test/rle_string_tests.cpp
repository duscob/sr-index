//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>

#include "sr-index/rle_string.hpp"

using BWT = std::string;
using Runs = std::vector<sri::StringRun>;

namespace sri {
void PrintTo(const StringRun &run, std::ostream *os) {
  *os << "{run: " << run.run << "; char: " << run.c << "; range: [" << run.range.first << "; " << run.range.second
      << "]}";
}
}

class RLEString_Tests : public testing::TestWithParam<std::tuple<BWT>> {};

TEST_P(RLEString_Tests, access_original) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

TEST_P(RLEString_Tests, access_byte_new) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<sdsl::wt_huff<>> bwt_rle(bwt.begin(), bwt.end());

  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

TEST_P(RLEString_Tests, access_int_new) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<sdsl::wt_huff_int<>> bwt_rle(bwt.begin(), bwt.end());

  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RLEString_Tests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}),
        std::make_tuple(BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'}),
        std::make_tuple(BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 'a', 'b', 'd', 'b', 'd', 'c', 'd'})
    )
);

class RLEStringRuns_Tests : public testing::TestWithParam<std::tuple<BWT, sri::range_t, Runs>> {};

TEST_P(RLEStringRuns_Tests, break_in_runs) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &range = std::get<1>(GetParam());
  auto runs_in_range = bwt_rle.break_in_runs(range);

  const auto &e_runs_in_range = std::get<2>(GetParam());
  EXPECT_THAT(runs_in_range, testing::ElementsAreArray(e_runs_in_range));
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RLEStringRuns_Tests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{0, 11},
                        Runs{sri::StringRun{0, 4, {0, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{4, 11},
                        Runs{sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{1, 11},
                        Runs{sri::StringRun{0, 4, {1, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{1, 7},
                        Runs{sri::StringRun{0, 4, {1, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 8},
                        Runs{sri::StringRun{4, 2, {5, 8}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{6, 8},
                        Runs{sri::StringRun{4, 2, {6, 8}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 7},
                        Runs{sri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{6, 7},
                        Runs{sri::StringRun{4, 2, {6, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{7, 7},
                        Runs{sri::StringRun{4, 2, {7, 7}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{3, 3},
                        Runs{sri::StringRun{2, 4, {3, 3}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 10},
                        Runs{sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 10}}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{7, 9},
                        Runs{sri::StringRun{4, 2, {7, 8}}, sri::StringRun{5, 3, {9, 9}}})
    )
);
