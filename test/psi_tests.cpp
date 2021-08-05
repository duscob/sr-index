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

  auto psi = sri::ConstructPsi(bwt, alphabet_);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, store) {
  const auto &bwt = std::get<0>(GetParam());

  sri::ConstructPsi(bwt, alphabet_, config_);

  sdsl::int_vector<> psi;
  sdsl::load_from_cache(psi, sdsl::conf::KEY_PSI, config_);

  const auto &e_psi = std::get<1>(GetParam());
  EXPECT_THAT(psi, testing::ElementsAreArray(e_psi));
}

TEST_P(PsiTests, partial_psi) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore(alphabet_.C, e_psi);

  auto psi_select = sri::RandomAccessForCRefContainer(std::cref(psi_core.select_partial_psi));
  auto get_c = [this](auto tt_index) { return sri::GetCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, partial_psi_serialize) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core";

  {
    auto tmp_psi_core = sri::PsiCore(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCore<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_select = sri::RandomAccessForCRefContainer(std::cref(psi_core.select_partial_psi));
  auto get_c = [this](auto tt_index) { return sri::GetCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, partial_psi_sd_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore<sdsl::sd_vector<>>(alphabet_.C, e_psi);

  auto psi_select = sri::RandomAccessForCRefContainer(std::cref(psi_core.select_partial_psi));
  auto get_c = [this](auto tt_index) { return sri::GetCForSAIndex(this->alphabet_.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  for (int i = 0; i < e_psi.size(); ++i) {
    EXPECT_EQ(psi(i), e_psi[i]) << "psi failed at index " << i;
  }
}

TEST_P(PsiTests, partial_psi_rrr_vector) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore<sdsl::rrr_vector<>>(alphabet_.C, e_psi);

  auto psi_select = sri::RandomAccessForCRefContainer(std::cref(psi_core.select_partial_psi));
  auto get_c = [this](auto tt_index) { return sri::GetCForSAIndex(this->alphabet_.C, tt_index); };
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

TEST_P(LFOnPsiTests, lf) {
  const auto &e_psi = std::get<1>(GetParam());

  auto psi_core = sri::PsiCore(alphabet_.C, e_psi);

  auto psi_rank = sri::RandomAccessForCRefContainer(std::cref(psi_core.rank_partial_psi));
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet_.C));

  auto lf = sri::LFOnPsi(psi_rank, cumulative);

  const auto &range = std::get<2>(GetParam());
  Char c = std::get<3>(GetParam());

  auto new_range = lf(range, alphabet_.char2comp[c]);

  auto e_range = std::get<4>(GetParam());
  EXPECT_EQ(new_range, e_range);
}

TEST_P(LFOnPsiTests, lf_serialized_phi) {
  const auto &e_psi = std::get<1>(GetParam());
  auto key = "psi_core";

  {
    auto tmp_psi_core = sri::PsiCore(alphabet_.C, e_psi);
    sdsl::store_to_cache(tmp_psi_core, key, config_);
  }

  sri::PsiCore<> psi_core;
  sdsl::load_from_cache(psi_core, key, config_);

  auto psi_rank = sri::RandomAccessForCRefContainer(std::cref(psi_core.rank_partial_psi));
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
