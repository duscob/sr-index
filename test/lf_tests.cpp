//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/21/21.
//

#include <gtest/gtest.h>

#include "sr-index/psi.h"
#include "sr-index/lf.h"
#include "sr-index/tools.h"

#include "base_tests.h"

class LFOnPsiTests : public BaseAlphabetTests, public testing::WithParamInterface<std::tuple<BWT, Psi, Range, Char, Range>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(LFOnPsiTests, psi_core_bv) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, psi_core_bv_serialized) {
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

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, psi_core_rle_serialized) {
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

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, psi_core_bv_is_trivial) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto is_trivial = psi_core.exist(alphabet_.char2comp[c], range.first);

  EXPECT_EQ(is_trivial, bwt_buf_[range.first] == c);
}

TEST_P(LFOnPsiTests, psi_core_rle_is_trivial) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto is_trivial = psi_core.exist(alphabet_.char2comp[c], range.first);

  EXPECT_EQ(is_trivial, bwt_buf_[range.first] == c);
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    LFOnPsiTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 12},
                        1,
                        Range{0, 1}), // 1 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 12},
                        2,
                        Range{1, 5}), // 2 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 12},
                        3,
                        Range{5, 9}), // 3 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 12},
                        4,
                        Range{9, 12}), // 4 before empty

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 5},
                        1,
                        Range{0, 1}), // 1 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 5},
                        2,
                        Range{1, 1}), // 2 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 5},
                        3,
                        Range{5, 6}), // 3 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 5},
                        4,
                        Range{10, 12}), // 4 before 2

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{10, 12},
                        2,
                        Range{5, 5}), // 2 before 4
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{10, 12},
                        3,
                        Range{7, 9}), // 3 before 4

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 9},
                        3,
                        Range{5, 6}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 10},
                        3,
                        Range{5, 7}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 11},
                        3,
                        Range{5, 8}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{2, 5},
                        3,
                        Range{5, 6}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{2, 9},
                        3,
                        Range{5, 6}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{2, 10},
                        3,
                        Range{5, 7}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{3, 5},
                        3,
                        Range{6, 6}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{3, 10},
                        3,
                        Range{6, 7}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{9, 12},
                        3,
                        Range{6, 9})
    )
);
