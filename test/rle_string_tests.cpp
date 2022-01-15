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

class AccessTests : public testing::TestWithParam<std::tuple<BWT>> {};

TEST_P(AccessTests, rle_string) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  EXPECT_EQ(bwt_rle.size(), bwt.size());
  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

TEST_P(AccessTests, RLEString_byte) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<sdsl::wt_huff<>> bwt_rle(bwt.begin(), bwt.end());

  EXPECT_EQ(bwt_rle.size(), bwt.size());
  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

TEST_P(AccessTests, RLEString_int) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<sdsl::wt_huff_int<>> bwt_rle(bwt.begin(), bwt.end());

  EXPECT_EQ(bwt_rle.size(), bwt.size());
  for (int i = 0; i < bwt.size(); ++i) {
    EXPECT_EQ(bwt_rle[i], bwt[i]) << "Failed at " << i;
  }
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    AccessTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}),
        std::make_tuple(BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'}),
        std::make_tuple(BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 'a', 'b', 'd', 'b', 'd', 'c', 'd'})
    )
);

using Char = unsigned char;

class SelectTests : public testing::TestWithParam<std::tuple<BWT, std::size_t, Char, std::size_t>> {};

TEST_P(SelectTests, rle_string) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &rnk = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto pos = bwt_rle.select(rnk - 1, c);

  const auto &e_pos = std::get<3>(GetParam());
  EXPECT_EQ(pos, e_pos);
}

TEST_P(SelectTests, StringRLE) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<> bwt_rle(bwt.begin(), bwt.end());

  const auto &rnk = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto pos = bwt_rle.select(rnk, c);

  const auto &e_pos = std::get<3>(GetParam());
  EXPECT_EQ(pos, e_pos);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    SelectTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 1, 4),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 2, 5),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 2, 6),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 2, 8),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 3, 2),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 3, 9),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 3, 11),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 4, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 4, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 4, 3)
    )
);

class RankTests : public testing::TestWithParam<std::tuple<BWT, std::size_t, Char, std::size_t>> {};

TEST_P(RankTests, rle_string) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &pos = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto rnk = bwt_rle.rank(pos, c);

  const auto &e_rnk = std::get<3>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
}

TEST_P(RankTests, StringRLE) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<> bwt_rle(bwt.begin(), bwt.end());

  const auto &pos = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto rnk = bwt_rle.rank(pos, c);

  const auto &e_rnk = std::get<3>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RankTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 1, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 1, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 5, 1, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 1, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 2, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 5, 2, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 6, 2, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 8, 2, 3),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 9, 2, 4),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 2, 4),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 3, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 3, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 9, 3, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 3, 4),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 0, 4, 0),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 4, 1),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 4, 2),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 4, 2),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 4, 3),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 4, 3)
    )
);

class SplitInRunsTests : public testing::TestWithParam<std::tuple<BWT, sri::range_t, Runs>> {};

TEST_P(SplitInRunsTests, rle_string) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &range = std::get<1>(GetParam());
  auto runs_in_range = bwt_rle.break_in_runs(range);

  const auto &e_runs_in_range = std::get<2>(GetParam());
  EXPECT_THAT(runs_in_range, testing::ElementsAreArray(e_runs_in_range));
}

TEST_P(SplitInRunsTests, RLEString) {
  const auto &bwt = std::get<0>(GetParam());
  sri::StringRLE<> bwt_rle(bwt.begin(), bwt.end());

  const auto &range = std::get<1>(GetParam());
  std::vector<sri::StringRun> runs_in_range;
  auto report = [&runs_in_range](auto tt_idx, auto tt_c, auto tt_start, auto tt_end) {
    runs_in_range.emplace_back(sri::StringRun{tt_idx, tt_c, sri::range_t{tt_start, tt_end - 1}});
  };
  bwt_rle.splitInRuns(range.first, range.second + 1, report);

  const auto &e_runs_in_range = std::get<2>(GetParam());
  EXPECT_THAT(runs_in_range, testing::ElementsAreArray(e_runs_in_range));
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    SplitInRunsTests,
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
