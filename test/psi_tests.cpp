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

using BWT = sdsl::int_vector<8>;

class BasePsiTests : public testing::Test {
 protected:
  void SetUp(const BWT &bwt) {
    n_ = bwt.size();

    sdsl::store_to_cache(bwt, key_tmp_bwt_, config_);

    bwt_buf_ = sdsl::int_vector_buffer<8>(cache_file_name(key_tmp_bwt_, config_));

    alphabet_ = sdsl::byte_alphabet(bwt_buf_, n_);
  }

  void TearDown() override {
    sdsl::util::delete_all_files(config_.file_map);
  }

  sdsl::cache_config config_;
  std::string key_tmp_bwt_ = sdsl::util::to_string(sdsl::util::pid()) + "_" + sdsl::util::to_string(sdsl::util::id());

  std::size_t n_ = 0;
  sdsl::int_vector_buffer<8> bwt_buf_;
  sdsl::byte_alphabet alphabet_;
};

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

using Psi = sdsl::int_vector<>;

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

TEST_P(PsiTests, psi_core) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core";

  {
    auto tmp_psi_core = sri::PsiCore(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCore<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_sd_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore<sdsl::sd_vector<>>(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, psi_core_rrr_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore<sdsl::rrr_vector<>>(alphabet_.C, e_psi);

  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };
  auto get_c = [this](auto tt_index) { return sri::computeCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
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

using Range = std::pair<std::size_t, std::size_t>;
using Char = unsigned char;

class LFOnPsiTests : public BasePsiTests, public testing::WithParamInterface<std::tuple<BWT, Psi, Range, Char, Range>> {
 protected:
  void SetUp() override {
    const auto &bwt = std::get<0>(GetParam());
    BasePsiTests::SetUp(bwt);
  }
};

TEST_P(LFOnPsiTests, psi_core) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore(alphabet_.C, e_psi);

  auto psi_rank = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.rank(tt_c, tt_rnk); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, psi_core_serialized) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core";

  {
    auto tmp_psi_core = sri::PsiCore(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCore<> psi_core;
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

INSTANTIATE_TEST_SUITE_P(
    Psi,
    LFOnPsiTests,
    testing::Values(
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 11},
                        1,
                        Range{0, 0}), // 1 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 11},
                        2,
                        Range{1, 4}), // 2 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 11},
                        3,
                        Range{5, 8}), // 3 before empty
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{0, 11},
                        4,
                        Range{9, 11}), // 4 before empty

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 4},
                        1,
                        Range{0, 0}), // 1 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 4},
                        2,
                        Range{1, 0}), // 2 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 4},
                        3,
                        Range{5, 5}), // 3 before 2
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{1, 4},
                        4,
                        Range{10, 11}), // 4 before 2

        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{10, 11},
                        2,
                        Range{1, 0}), // 2 before 42
        std::make_tuple(BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        Range{10, 11},
                        3,
                        Range{7, 8}) // 3 before 42
    )
);
