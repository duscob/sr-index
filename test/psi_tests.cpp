//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/5/21.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>

#include "sr-index/psi.h"
#include "sr-index/tools.h"

#include "psi_base_tests.h"

using Cumulative = sdsl::int_vector<64>;

class AlphabetTests : public BasePsiTests, public testing::WithParamInterface<std::tuple<BWT, Cumulative>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BasePsiTests::SetUp(bwt);
  }
};

TEST_P(AlphabetTests, construct) {
  const auto &cumulative = this->alphabet_.C;

  const auto &e_cumulative = std::get<1>(GetParam());
  EXPECT_THAT(cumulative, testing::ElementsAreArray(e_cumulative));
}

INSTANTIATE_TEST_SUITE_P(
    Alphabet,
    AlphabetTests,
    testing::Values(
        std::make_tuple(
            BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
            Cumulative{0, 1, 5, 9, 12}),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
            Cumulative{0, 1, 5, 9, 12}),
        std::make_tuple(
            BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 'a', 'b', 'd', 'b', 'd', 'c', 'd'},
            Cumulative{0, 1, 5, 6, 12, 16})
    )
);

class PsiTests : public BasePsiTests, public testing::WithParamInterface<std::tuple<BWT, Psi>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BasePsiTests::SetUp(bwt);
  }
};

TEST_P(PsiTests, construct) {
  const auto &bwt = std::get<0>(GetParam());

  auto psi = sri::constructPsi(bwt, alphabet_);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, store) {
  const auto &bwt = std::get<0>(GetParam());

  sri::constructPsi(bwt, alphabet_, config_);

  sdsl::int_vector<> psi;
  sdsl::load_from_cache(psi, sdsl::conf::KEY_PSI, config_);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, psi_core_bv) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  EXPECT_EQ(psi_core.getFirstBWTSymbol(), alphabet_.char2comp[bwt_buf_[0]]);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_bv_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core_bv";

  {
    auto tmp_psi_core = sri::PsiCoreBV(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCoreBV<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_bv_sd_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV<sdsl::sd_vector<>>(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_bv_rrr_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreBV<sdsl::rrr_vector<>>(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  EXPECT_EQ(psi_core.getFirstBWTSymbol(), alphabet_.char2comp[bwt_buf_[0]]);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_rle_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core_rle";

  {
    auto tmp_psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCoreRLE<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_rle_rank_run) {
  const auto &bwt = std::get<0>(GetParam());
  sdsl::bit_vector bwt_run_start(bwt.size(), 0);
  bwt_run_start[0] = true;
  for (int i = 1; i < bwt.size(); ++i) {
    if (bwt[i - 1] != bwt[i]) {
      bwt_run_start[i] = true;
    }
  }
  sdsl::bit_vector::rank_1_type bwt_run_start_rank(&bwt_run_start);

  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  for (int i = 0; i <= e_psi.size(); ++i) {
    EXPECT_EQ(psi_core.rankRun(i), bwt_run_start_rank(i)) << "psi failed at index " << i;
  }
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    PsiTests,
    testing::Values(
        std::make_tuple(
            BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}),
        std::make_tuple(
            BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 'a', 'b', 'd', 'b', 'd', 'c', 'd'},
            Psi{9, 7, 8, 10, 12, 14, 0, 5, 6, 11, 13, 15, 1, 2, 3, 4})
    )
);

using Ranges = std::vector<Range>;

class SplitInRunsWithPsiTests
    : public BasePsiTests, public testing::WithParamInterface<std::tuple<BWT, Psi, Range, Ranges>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BasePsiTests::SetUp(bwt);
  }
};

TEST_P(SplitInRunsWithPsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &range = std::get<2>(GetParam());

  auto ranges = psi_core.splitInRuns(range.first, range.second);

  auto e_ranges = std::get<3>(GetParam());
  EXPECT_THAT(ranges, testing::ElementsAreArray(e_ranges));
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    SplitInRunsWithPsiTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 12},
                        Ranges{{0, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 9}, {9, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 5},
                        Ranges{{1, 2}, {2, 3}, {3, 4}, {4, 5}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{2, 5},
                        Ranges{{2, 3}, {3, 4}, {4, 5}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{3, 5},
                        Ranges{{3, 4}, {4, 5}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{4, 5},
                        Ranges{{4, 5}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{5, 9},
                        Ranges{{5, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{6, 9},
                        Ranges{{6, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{7, 9},
                        Ranges{{7, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{8, 9},
                        Ranges{{8, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{9, 12},
                        Ranges{{9, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{10, 12},
                        Ranges{{10, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{9, 11},
                        Ranges{{9, 11}})
    )
);

class ComputeForwardRunsWithPsiTests
    : public BasePsiTests, public testing::WithParamInterface<std::tuple<BWT, Psi, Char, Range, Ranges>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BasePsiTests::SetUp(bwt);
  }
};

TEST_P(ComputeForwardRunsWithPsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &c = std::get<2>(GetParam());
  const auto &range = std::get<3>(GetParam());

  auto ranges = psi_core.computeForwardRuns(c, range.first, range.second);

  auto e_ranges = std::get<4>(GetParam());
  EXPECT_THAT(ranges, testing::ElementsAreArray(e_ranges));
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    ComputeForwardRunsWithPsiTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        0,
                        Range{1, 2},
                        Ranges{{4, 5}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        1,
                        Range{1, 5},
                        Ranges{{5, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        2,
                        Range{1, 5},
                        Ranges{{2, 3}, {9, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        3,
                        Range{1, 4},
                        Ranges{{0, 2}, {3, 4}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        1,
                        Range{1, 2},
                        Ranges{{5, 6}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        1,
                        Range{2, 5},
                        Ranges{{6, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        1,
                        Range{3, 5},
                        Ranges{{7, 9}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        1,
                        Range{3, 4},
                        Ranges{{7, 8}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        2,
                        Range{2, 5},
                        Ranges{{9, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        2,
                        Range{3, 5},
                        Ranges{{10, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        2,
                        Range{4, 5},
                        Ranges{{11, 12}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        2,
                        Range{3, 4},
                        Ranges{{10, 11}}),
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        3,
                        Range{2, 4},
                        Ranges{{1, 2}, {3, 4}})
    )
);
