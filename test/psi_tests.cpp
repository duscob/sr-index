//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/5/21.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>

#include "sr-index/psi.h"
#include "sr-index/rle_string.hpp"
#include "sr-index/tools.h"

#include "base_tests.h"

using Cumulative = sdsl::int_vector<64>;

class AlphabetTests : public BaseAlphabetTests, public testing::WithParamInterface<std::tuple<BWT, Cumulative>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
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
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
            Cumulative{0, 1, 5, 9, 12}),
        std::make_tuple(
            BWT{4, 4, 3, 4, 0, 2, 2, 2, 2, 3, 3, 3},
            Cumulative{0, 1, 5, 9, 12}),
        std::make_tuple(
            BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 0, 'b', 'd', 'b', 'd', 'c', 'd'},
            Cumulative{0, 1, 5, 6, 12, 16})
    )
);

class PsiTests : public BaseAlphabetTests, public testing::WithParamInterface<std::tuple<BWT, Psi>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(PsiTests, construct) {
  const auto &bwt = std::get<0>(GetParam());
  auto get_bwt_symbol = [&bwt, this](auto tt_i) { return this->alphabet_.char2comp[bwt[tt_i]]; };

  auto psi = sri::constructPsi(get_bwt_symbol, alphabet_.C);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, construct_bwt_rle) {
  const auto &bwt = std::get<0>(GetParam());
  auto get_bwt_symbol = [&bwt, this](auto tt_i) { return this->alphabet_.char2comp[bwt[tt_i]]; };
  auto bwt_s = sdsl::random_access_container(get_bwt_symbol, bwt.size());
  sri::RLEString<> bwt_rle(bwt_s.begin(), bwt_s.end());
  auto get_bwt_rle_symbol = [&bwt_rle](auto tt_i) { return bwt_rle[tt_i]; };

  auto psi = sri::constructPsi(get_bwt_rle_symbol, alphabet_.C);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, store) {
  const auto &bwt = std::get<0>(GetParam());
  auto get_bwt_symbol = [&bwt, this](size_t tt_i) { return this->alphabet_.char2comp[bwt[tt_i]]; };

  sri::constructPsi(get_bwt_symbol, alphabet_.C, config_);

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
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}),
        std::make_tuple(
            BWT{4, 4, 3, 4, 0, 2, 2, 2, 2, 3, 3, 3},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}),
        std::make_tuple(
            BWT{'d', 'e', 'e', 'e', 'e', 'd', 'd', 'b', 'b', 0, 'b', 'd', 'b', 'd', 'c', 'd'},
            Psi{9, 7, 8, 10, 12, 14, 0, 5, 6, 11, 13, 15, 1, 2, 3, 4})
    )
);

using Value = std::size_t;
struct DataRank {
  std::size_t rank = 0;

  struct RunRank {
    std::size_t start = 0;
    std::size_t end = 0;
    std::size_t rank = 0;
  } run;

  bool operator==(const DataRank &rhs) const {
    return rank == rhs.rank && run.start == rhs.run.start && run.end == rhs.run.end && run.rank == rhs.run.rank;
  }
};

void PrintTo(const DataRank &t_item, std::ostream *t_os) {
  *t_os << "{" << t_item.rank << "; {" << t_item.run.start << "; " << t_item.run.end << "; " << t_item.run.rank << "}}";
}

class RankPsiTests : public BaseAlphabetTests,
                     public testing::WithParamInterface<std::tuple<BWT, Psi, std::tuple<Char, Value, DataRank>>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(RankPsiTests, psi_core_rle_rank_simple) {
  const auto &e_psi = std::get<1>(GetParam());
  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &c = std::get<0>(item);
  const auto &value = std::get<1>(item);

  auto rank = psi_core.rank(c, value);

  const auto &e_data_rank = std::get<2>(item);
  EXPECT_EQ(rank, e_data_rank.rank);
}

TEST_P(RankPsiTests, psi_core_rle_rank_extra) {
  const auto &e_psi = std::get<1>(GetParam());
  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &c = std::get<0>(item);
  const auto &value = std::get<1>(item);

  DataRank data_rank;
  auto report =
      [&data_rank](const auto &tt_rank, const auto &tt_run_start, const auto &tt_run_end, const auto &tt_run_rank) {
        data_rank = DataRank{tt_rank, {tt_run_start, tt_run_end, tt_run_rank}};
      };
  psi_core.rank(c, value, report);

  const auto &e_data_rank = std::get<2>(item);
  EXPECT_EQ(data_rank, e_data_rank);
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    RankPsiTests,
    testing::Combine(
        testing::Values(
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'}
        ),
        testing::Values(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}
        ),
        testing::Values(
            std::make_tuple(0, 0, DataRank{0, {4, 5, 0}}),
            std::make_tuple(0, 3, DataRank{0, {4, 5, 0}}),
            std::make_tuple(0, 4, DataRank{0, {4, 5, 0}}),
            std::make_tuple(0, 5, DataRank{1, {12, 12, 1}}),
            std::make_tuple(0, 6, DataRank{1, {12, 12, 1}}),
            std::make_tuple(0, 12, DataRank{1, {12, 12, 1}}),

            std::make_tuple(1, 0, DataRank{0, {5, 9, 0}}),
            std::make_tuple(1, 4, DataRank{0, {5, 9, 0}}),
            std::make_tuple(1, 5, DataRank{0, {5, 9, 0}}),
            std::make_tuple(1, 6, DataRank{1, {5, 9, 0}}),
            std::make_tuple(1, 7, DataRank{2, {5, 9, 0}}),
            std::make_tuple(1, 8, DataRank{3, {5, 9, 0}}),
            std::make_tuple(1, 9, DataRank{4, {12, 12, 1}}),
            std::make_tuple(1, 10, DataRank{4, {12, 12, 1}}),
            std::make_tuple(1, 12, DataRank{4, {12, 12, 1}}),

            std::make_tuple(2, 0, DataRank{0, {2, 3, 0}}),
            std::make_tuple(2, 1, DataRank{0, {2, 3, 0}}),
            std::make_tuple(2, 2, DataRank{0, {2, 3, 0}}),
            std::make_tuple(2, 3, DataRank{1, {9, 12, 1}}),
            std::make_tuple(2, 8, DataRank{1, {9, 12, 1}}),
            std::make_tuple(2, 9, DataRank{1, {9, 12, 1}}),
            std::make_tuple(2, 10, DataRank{2, {9, 12, 1}}),
            std::make_tuple(2, 11, DataRank{3, {9, 12, 1}}),
            std::make_tuple(2, 12, DataRank{4, {9, 12, 1}}),

            std::make_tuple(3, 0, DataRank{0, {0, 2, 0}}),
            std::make_tuple(3, 1, DataRank{1, {0, 2, 0}}),
            std::make_tuple(3, 2, DataRank{2, {3, 4, 1}}),
            std::make_tuple(3, 3, DataRank{2, {3, 4, 1}}),
            std::make_tuple(3, 4, DataRank{3, {12, 12, 2}}),
            std::make_tuple(3, 12, DataRank{3, {12, 12, 2}})
        )
    )
);

using Ranges = std::vector<Range>;

class SplitInRunsWithPsiTests : public BaseAlphabetTests,
                                public testing::WithParamInterface<std::tuple<BWT, Psi, std::tuple<Range, Ranges>>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(SplitInRunsWithPsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &range = std::get<0>(item);

  auto ranges = psi_core.splitInRuns(range.first, range.second);

  auto e_ranges = std::get<1>(item);
  EXPECT_THAT(ranges, testing::ElementsAreArray(e_ranges));
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    SplitInRunsWithPsiTests,
    testing::Combine(
        testing::Values(
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'}
        ),
        testing::Values(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}
        ),
        testing::Values(
            std::make_tuple(Range{0, 12}, Ranges{{0, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 9}, {9, 12}}),
            std::make_tuple(Range{1, 5}, Ranges{{1, 2}, {2, 3}, {3, 4}, {4, 5}}),
            std::make_tuple(Range{2, 5}, Ranges{{2, 3}, {3, 4}, {4, 5}}),
            std::make_tuple(Range{3, 5}, Ranges{{3, 4}, {4, 5}}),
            std::make_tuple(Range{4, 5}, Ranges{{4, 5}}),
            std::make_tuple(Range{5, 9}, Ranges{{5, 9}}),
            std::make_tuple(Range{6, 9}, Ranges{{6, 9}}),
            std::make_tuple(Range{7, 9}, Ranges{{7, 9}}),
            std::make_tuple(Range{8, 9}, Ranges{{8, 9}}),
            std::make_tuple(Range{9, 12}, Ranges{{9, 12}}),
            std::make_tuple(Range{10, 12}, Ranges{{10, 12}}),
            std::make_tuple(Range{9, 11}, Ranges{{9, 11}})
        )
    )
);

class ComputeForwardRunsWithPsiTests : public BaseAlphabetTests,
                                       public testing::WithParamInterface<
                                           std::tuple<BWT, Psi, std::tuple<Char, Range, Ranges>>
                                       > {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BaseAlphabetTests::SetUp(bwt);
  }
};

TEST_P(ComputeForwardRunsWithPsiTests, psi_core_rle) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCoreRLE(alphabet_.C, e_psi);

  const auto &item = std::get<2>(GetParam());
  const auto &c = std::get<0>(item);
  const auto &range = std::get<1>(item);

  auto ranges = psi_core.computeForwardRuns(c, range.first, range.second);

  auto e_ranges = std::get<2>(item);
  EXPECT_THAT(ranges, testing::ElementsAreArray(e_ranges));
}

INSTANTIATE_TEST_SUITE_P(
    Psi,
    ComputeForwardRunsWithPsiTests,
    testing::Combine(
        testing::Values(
            BWT{'c', 'c', 'b', 'c', 0, 'a', 'a', 'a', 'a', 'b', 'b', 'b'}
        ),
        testing::Values(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}
        ),
        testing::Values(
            std::make_tuple(0, Range{1, 2}, Ranges{{4, 5}}),
            std::make_tuple(1, Range{1, 5}, Ranges{{5, 9}}),
            std::make_tuple(2, Range{1, 5}, Ranges{{2, 3}, {9, 12}}),
            std::make_tuple(3, Range{1, 4}, Ranges{{0, 2}, {3, 4}}),
            std::make_tuple(1, Range{1, 2}, Ranges{{5, 6}}),
            std::make_tuple(1, Range{2, 5}, Ranges{{6, 9}}),
            std::make_tuple(1, Range{3, 5}, Ranges{{7, 9}}),
            std::make_tuple(1, Range{3, 4}, Ranges{{7, 8}}),
            std::make_tuple(2, Range{2, 5}, Ranges{{9, 12}}),
            std::make_tuple(2, Range{3, 5}, Ranges{{10, 12}}),
            std::make_tuple(2, Range{4, 5}, Ranges{{11, 12}}),
            std::make_tuple(2, Range{3, 4}, Ranges{{10, 11}}),
            std::make_tuple(3, Range{2, 4}, Ranges{{1, 2}, {3, 4}})
        )
    )
);
