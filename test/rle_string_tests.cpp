//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/26/20.
//

#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>

#include "sr-index/rle_string.hpp"

using String = std::string;
using Runs = std::vector<sri::StringRun>;

namespace sri {
void PrintTo(const StringRun &run, std::ostream *os) {
  *os << "{run: " << run.run << "; char: " << run.c << "; range: [" << run.range.first << "; " << run.range.second
      << "]}";
}
}

class AccessTests : public testing::TestWithParam<std::tuple<String>> {};

TEST_P(AccessTests, rle_string) {
  const auto &str = std::get<0>(GetParam());
  sri::rle_string<> rle_str(str);

  EXPECT_EQ(rle_str.size(), str.size());
  for (int i = 0; i < str.size(); ++i) {
    EXPECT_EQ(rle_str[i], str[i]) << "Failed at " << i;
  }
}

TEST_P(AccessTests, RLEString_byte) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<sdsl::wt_huff<>> rle_str(str.begin(), str.end());

  EXPECT_EQ(rle_str.size(), str.size());
  for (int i = 0; i < str.size(); ++i) {
    EXPECT_EQ(rle_str[i], str[i]) << "Failed at " << i;
  }
}

TEST_P(AccessTests, RLEString_int) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<sdsl::wt_huff_int<>> rle_str(str.begin(), str.end());

  EXPECT_EQ(rle_str.size(), str.size());
  for (int i = 0; i < str.size(); ++i) {
    EXPECT_EQ(rle_str[i], str[i]) << "Failed at " << i;
  }
}

TEST_P(AccessTests, RLEString_io) {
  const auto &str = std::get<0>(GetParam());
  sdsl::cache_config config;
  auto key = "rle_bwt";

  {
    sri::RLEString<> tmp_rle_str(str.begin(), str.end());
    sdsl::store_to_cache(tmp_rle_str, key, config);
  }

  sri::RLEString<> rle_str;
  sdsl::load_from_cache(rle_str, key, config);

  EXPECT_EQ(rle_str.size(), str.size());
  for (int i = 0; i < str.size(); ++i) {
    EXPECT_EQ(rle_str[i], str[i]) << "Failed at " << i;
  }

  sdsl::util::delete_all_files(config.file_map);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    AccessTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}),
        std::make_tuple(String{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'}),
        std::make_tuple(String{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 'a', 'b', 'd', 'b', 'd', 'c', 'd'})
    )
);

using Char = unsigned char;

class SelectTests : public testing::TestWithParam<std::tuple<String, std::size_t, Char, std::size_t>> {};

TEST_P(SelectTests, rle_string) {
  const auto &str = std::get<0>(GetParam());
  sri::rle_string<> rle_str(str);

  const auto &rnk = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto pos = rle_str.select(rnk - 1, c);

  const auto &e_pos = std::get<3>(GetParam());
  EXPECT_EQ(pos, e_pos);
}

TEST_P(SelectTests, RLEString) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &rnk = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto pos = rle_str.select(rnk, c);

  const auto &e_pos = std::get<3>(GetParam());
  EXPECT_EQ(pos, e_pos);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    SelectTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 1, 4),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 2, 5),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 2, 6),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 2, 8),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 3, 2),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 3, 9),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 3, 11),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 4, 1),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 4, 3)
    )
);

class RankTests : public testing::TestWithParam<std::tuple<String, std::size_t, Char, std::size_t, std::size_t, bool>> {
};

TEST_P(RankTests, rle_string) {
  const auto &str = std::get<0>(GetParam());
  sri::rle_string<> rle_str(str);

  const auto &pos = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto rnk = rle_str.rank(pos, c);

  const auto &e_rnk = std::get<3>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
}

TEST_P(RankTests, RLEString) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &pos = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto rnk = rle_str.rank(pos, c);

  const auto &e_rnk = std::get<3>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
}

TEST_P(RankTests, RLEString_report) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &pos = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  std::size_t rnk;
  std::size_t run_rnk;
  bool contained;
  auto report = [&rnk, &run_rnk, &contained](auto tt_rnk, auto tt_run_rnk, auto tt_contained) {
    rnk = tt_rnk;
    run_rnk = tt_run_rnk;
    contained = tt_contained;
  };
  rle_str.rank(pos, c, report);

  const auto &e_rnk = std::get<3>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
  const auto &e_run_rnk = std::get<4>(GetParam());
  EXPECT_EQ(run_rnk, e_run_rnk);
  const auto &e_contained = std::get<5>(GetParam());
  EXPECT_EQ(contained, e_contained);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RankTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 1, 0, 0, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 1, 0, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 5, 1, 1, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 1, 1, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 2, 0, 0, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 5, 2, 0, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 6, 2, 1, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 8, 2, 3, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 9, 2, 4, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 2, 4, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 3, 0, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 3, 1, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 9, 3, 1, 2, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 3, 4, 2, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 0, 4, 0, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 4, 1, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 4, 2, 1, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 4, 2, 2, true),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 4, 3, 2, false),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 12, 4, 3, 2, false),

        std::make_tuple(String{4, 4, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3}, 4, 5, 0, 1, true),
        std::make_tuple(String{4, 4, 3, 4, 5, 2, 2, 2, 2, 3, 3, 3}, 4, 1, 0, 0, false)
    )
);

class RankRawTests : public testing::TestWithParam<
    std::tuple<String, std::size_t, std::size_t, Char, std::size_t, std::size_t>> {
};

TEST_P(RankRawTests, RLEString_report) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &pos = std::get<1>(GetParam());
  std::size_t rnk;
  uint8_t c;
  std::size_t run_rnk;
  std::size_t symbol_run_rnk;
  auto report = [&rnk, &c, &run_rnk, &symbol_run_rnk](auto tt_rnk, auto tt_c, auto tt_run_rnk, auto tt_symbol_run_rnk) {
    rnk = tt_rnk;
    c = tt_c;
    run_rnk = tt_run_rnk;
    symbol_run_rnk = tt_symbol_run_rnk;
  };
  rle_str.rank(pos, report);

  const auto &e_rnk = std::get<2>(GetParam());
  EXPECT_EQ(rnk, e_rnk);
  const auto &e_c = std::get<3>(GetParam());
  EXPECT_EQ(c, e_c);
  const auto &e_run_rnk = std::get<4>(GetParam());
  EXPECT_EQ(run_rnk, e_run_rnk);
  const auto &e_symbol_run_rnk = std::get<5>(GetParam());
  EXPECT_EQ(symbol_run_rnk, e_symbol_run_rnk);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    RankRawTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 0, 0, 4, 0, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 1, 4, 0, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 0, 3, 1, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 3, 2, 4, 2, 1),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 4, 0, 1, 3, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 5, 0, 2, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 6, 1, 2, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 7, 2, 2, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 8, 3, 2, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 9, 1, 3, 5, 1),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 10, 2, 3, 5, 1),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 11, 3, 3, 5, 1)
    )
);

class SelectRunTests : public testing::TestWithParam<std::tuple<String, std::size_t, Char, std::size_t>> {};

TEST_P(SelectRunTests, RLEString) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &rnk = std::get<1>(GetParam());
  const auto &c = std::get<2>(GetParam());
  auto pos = rle_str.selectOnRuns(rnk, c);

  const auto &e_pos = std::get<3>(GetParam());
  EXPECT_EQ(pos, e_pos);
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    SelectRunTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 1, 3),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 2, 4),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 3, 1),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 3, 5),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 1, 4, 0),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, 2, 4, 2)
    )
);

class SplitInRunsTests : public testing::TestWithParam<std::tuple<String, sri::range_t, Runs>> {};

TEST_P(SplitInRunsTests, rle_string) {
  const auto &str = std::get<0>(GetParam());
  sri::rle_string<> rle_str(str);

  const auto &range = std::get<1>(GetParam());
  auto runs = rle_str.break_in_runs(range);

  const auto &e_runs = std::get<2>(GetParam());
  EXPECT_THAT(runs, testing::ElementsAreArray(e_runs));
}

TEST_P(SplitInRunsTests, RLEString) {
  const auto &str = std::get<0>(GetParam());
  sri::RLEString<> rle_str(str.begin(), str.end());

  const auto &range = std::get<1>(GetParam());
  std::vector<sri::StringRun> runs;
  auto report = [&runs, &range](auto tt_idx, auto tt_c, auto tt_start, auto tt_end) {
    runs.emplace_back(sri::StringRun{
        tt_idx, tt_c, sri::range_t{std::max(tt_start, range.first), std::min(range.second, tt_end - 1)}});
  };
  rle_str.splitInRuns(range.first, range.second + 1, report);

  const auto &e_runs = std::get<2>(GetParam());
  EXPECT_THAT(runs, testing::ElementsAreArray(e_runs));
}

INSTANTIATE_TEST_SUITE_P(
    RLEString,
    SplitInRunsTests,
    testing::Values(
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{0, 11},
                        Runs{sri::StringRun{0, 4, {0, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{4, 11},
                        Runs{sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{1, 11},
                        Runs{sri::StringRun{0, 4, {1, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 11}}}
        ),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{1, 7},
                        Runs{sri::StringRun{0, 4, {1, 1}}, sri::StringRun{1, 3, {2, 2}}, sri::StringRun{2, 4, {3, 3}},
                             sri::StringRun{3, 1, {4, 4}}, sri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 8},
                        Runs{sri::StringRun{4, 2, {5, 8}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{6, 8},
                        Runs{sri::StringRun{4, 2, {6, 8}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 7},
                        Runs{sri::StringRun{4, 2, {5, 7}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{6, 7},
                        Runs{sri::StringRun{4, 2, {6, 7}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{7, 7},
                        Runs{sri::StringRun{4, 2, {7, 7}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{3, 3},
                        Runs{sri::StringRun{2, 4, {3, 3}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{5, 10},
                        Runs{sri::StringRun{4, 2, {5, 8}}, sri::StringRun{5, 3, {9, 10}}}),
        std::make_tuple(String{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        sri::range_t{7, 9},
                        Runs{sri::StringRun{4, 2, {7, 8}}, sri::StringRun{5, 3, {9, 9}}})
    )
);
