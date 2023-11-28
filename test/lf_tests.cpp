//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/21/21.
//

#include <gtest/gtest.h>

#include "sr-index/rle_string.hpp"
#include "sr-index/psi.h"
#include "sr-index/lf.h"
#include "sr-index/tools.h"

#include "base_tests.h"

class LFTests : public BaseAlphabetTests,
                public testing::WithParamInterface<std::tuple<BWT, Psi, std::tuple<Range, Char, Range>>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(LFTests, rle_string) {
  const auto &e_psi = std::get<1>(GetParam());

  auto get_symbol = [this](auto tt_i) { return this->alphabet_.char2comp[this->bwt_buf_[tt_i]]; };
  auto bwt_s = sdsl::random_access_container(get_symbol, this->bwt_buf_.size());
  sri::RLEString<> bwt_rle(bwt_s.begin(), bwt_s.end());

  auto bwt_rank = [&bwt_rle](auto tt_c, auto tt_rnk) { return bwt_rle.rank(tt_rnk, tt_c); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(bwt_rank, cumulative);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<2>(item);
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFTests, psi_core_bv) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<2>(item);
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFTests, psi_core_bv_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core_bv";

  {
    auto tmp_psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCoreBV<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<2>(item);
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<2>(item);
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFTests, psi_core_rle_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core_rle";

  {
    auto tmp_psi_rle = sri::PsiCoreRLE(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_rle, key, config_);
  }

  sri::PsiCoreRLE<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<2>(item);
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFTests, psi_core_bv_is_trivial) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto is_trivial = psi_core.exist(alphabet_.char2comp[c], range.first);

  EXPECT_EQ(is_trivial, bwt_buf_[range.first] == c);
}

TEST_P(LFTests, psi_core_rle_is_trivial) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);
  Char c = std::get<1>(item);

  auto is_trivial = psi_core.exist(alphabet_.char2comp[c], range.first);

  EXPECT_EQ(is_trivial, bwt_buf_[range.first] == c);
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    LFTests,
    testing::Combine(
        testing::Values(
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'}
        ),
        testing::Values(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}
        ),
        testing::Values(
            std::make_tuple(Range{0, 12}, 0, Range{0, 1}), // 1 before empty
            std::make_tuple(Range{0, 12}, 'a', Range{1, 5}), // 2 before empty
            std::make_tuple(Range{0, 12}, 'b', Range{5, 9}), // 3 before empty
            std::make_tuple(Range{0, 12}, 'c', Range{9, 12}), // 4 before empty

            std::make_tuple(Range{1, 5}, 0, Range{0, 1}), // 1 before 2
            std::make_tuple(Range{1, 5}, 'a', Range{1, 1}), // 2 before 2
            std::make_tuple(Range{1, 5}, 'b', Range{5, 6}), // 3 before 2
            std::make_tuple(Range{1, 5}, 'c', Range{10, 12}), // 4 before 2

            std::make_tuple(Range{10, 12}, 'a', Range{5, 5}), // 2 before 4
            std::make_tuple(Range{10, 12}, 'b', Range{7, 9}), // 3 before 4

            std::make_tuple(Range{1, 9}, 'b', Range{5, 6}),
            std::make_tuple(Range{1, 10}, 'b', Range{5, 7}),
            std::make_tuple(Range{1, 11}, 'b', Range{5, 8}),
            std::make_tuple(Range{2, 5}, 'b', Range{5, 6}),
            std::make_tuple(Range{2, 9}, 'b', Range{5, 6}),
            std::make_tuple(Range{2, 10}, 'b', Range{5, 7}),
            std::make_tuple(Range{3, 5}, 'b', Range{6, 6}),
            std::make_tuple(Range{3, 10}, 'b', Range{6, 7}),
            std::make_tuple(Range{9, 12}, 'b', Range{6, 9})
        )
    )
);
