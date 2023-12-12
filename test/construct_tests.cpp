//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/2023.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/config.h"
#include "sr-index/r_csa.h"
#include "sr-index/sr_csa_psi.h"

#include "base_tests.h"

using PsiRunHead = IntVector;
using PsiRunHeadAsc = IntVector;
using PsiRunHeadInd = IntVector;
using PsiRunTail = IntVector;
using PsiRunTailAsc = IntVector;
using PsiRunTailAscLink = IntVector;
using BitVector = sdsl::bit_vector;
using Marks = BitVector;
using SampleRate = std::size_t;

class BaseConstructTests : public BaseConfigTests {
public:
  void compare(const std::string& tt_key, const sdsl::int_vector<>& tt_e_values) const {
    sdsl::int_vector<> values;
    load_from_cache(values, tt_key, config_);

    EXPECT_THAT(values, testing::ElementsAreArray(tt_e_values)) << "Key = " << tt_key;
  }

  void compare(const std::string& tt_key, const sdsl::sd_vector<>& tt_e_values) const {
    sdsl::sd_vector<> values;
    load_from_cache(values, tt_key, config_, true);

    EXPECT_THAT(values, testing::ElementsAreArray(tt_e_values)) << "Key = " << tt_key;
  }
};

class ConstructRCSATests : public BaseConstructTests,
                           public testing::WithParamInterface<std::tuple<
                             String, Psi, PsiRunHead, PsiRunTail, PsiRunTailAsc, PsiRunTailAscLink, Marks
                           >> {
protected:
  void SetUp() override {
    const auto& data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(ConstructRCSATests, construct) {
  using namespace sri::conf;
  sri::constructRCSAWithPsiRuns<8, sdsl::sd_vector<>>(config_.file_map[key_tmp_input_], config_);

  compare(config_.keys[kPsi][kBase], std::get<1>(GetParam()));
  compare(config_.keys[kPsi][kHead][kTextPos], std::get<2>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPos], std::get<3>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPosAsc][kIdx], std::get<4>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPosAsc][kLink], std::get<5>(GetParam()));

  compare(config_.keys[kPsi][kTail][kTextPos], std::get<6>(GetParam()));
}

INSTANTIATE_TEST_SUITE_P(
  Basic,
  ConstructRCSATests,
  testing::Values(
    std::make_tuple(
      String{"alabaralaalabarda"},
      Psi{6, 0, 7, 10, 11, 13, 14, 15, 16, 17, 8, 9, 1, 2, 3, 4, 5, 12},
      // Psi{6,    0,    7,    10,   11,   13,   14,   15,   16,   17,   8,    9,    1,    2,    3,    4,    5,    12},
      // SA{17,    16,   8,    2,    11,   6,    0,    9,    4,    13,   3,    12,   15,   7,    1,    10,   5,    14}
      // BWT{'a',  'd',  'l',  'l',  'l',  'r',  '$',  'a',  'b',  'b',  'a',  'a',  'r',  'a',  'a',  'a',  'a',  'a'},
      PsiRunHead{17, 16, 8, 2, 6, 3, 15, 7, 5, 14},
      PsiRunTail{17, 16, 8, 11, 13, 12, 15, 10, 5, 14},
      PsiRunTailAsc{8, 2, 7, 3, 5, 4, 9, 6, 1, 0},
      PsiRunTailAscLink{9, 3, 8, 4, 6, 5, 0, 7, 2, 1},
      Marks{0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}
    ),
    std::make_tuple(
      String{"abcabcababc"},
      Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
      // Psi{4,    5,    6,    7,    8,    2,    9,    10,   11,   0,    1,    3},
      // SA{11,    6,    8,    3,    0,    7,    9,    4,    1,    10,   5,    2},
      // BWT{'c',  'c',  'b',  'c',  '$',  'a',  'a',  'a',  'a',  'b',  'b',  'b'},
      PsiRunHead{11, 6, 7, 9, 10, 2},
      PsiRunTail{11, 0, 7, 1, 5, 2},
      PsiRunTailAsc{1, 3, 5, 4, 2, 0},
      PsiRunTailAscLink{2, 4, 0, 5, 3, 1},
      Marks{1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1}
    )
  )
);

class ConstructSRCSATests : public BaseConstructTests,
                            public testing::WithParamInterface<std::tuple<
                              String, SampleRate, PsiRunHeadAsc, PsiRunHeadInd, PsiRunHead>> {
protected:
  void SetUp() override {
    const auto& data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(ConstructSRCSATests, construct) {
  using namespace sri::conf;

  sri::constructRCSAWithPsiRuns<8, sdsl::sd_vector<>>(config_.file_map[key_tmp_input_], config_);
  const auto& subsample_rate = std::get<1>(GetParam());
  sri::constructSrCSACommonsWithPsiRuns<8, sdsl::sd_vector<>>(subsample_rate, config_);

  const auto prefix = std::to_string(subsample_rate) + "_";

  compare(config_.keys[kPsi][kHead][kTextPosAsc][kIdx], std::get<2>(GetParam()));
  compare(prefix + to_string(config_.keys[kPsi][kHead][kIdx]), std::get<3>(GetParam()));
  compare(prefix + to_string(config_.keys[kPsi][kHead][kTextPos]), std::get<4>(GetParam()));
}

INSTANTIATE_TEST_SUITE_P(
  Basic,
  ConstructSRCSATests,
  testing::Values(
    std::make_tuple(
      String{"alabaralaalabarda"},
      4,
      PsiRunHeadAsc{3, 5, 8, 4, 7, 2, 9, 6, 1, 0},
      PsiRunHeadInd{0, 2, 3, 4, 9},
      PsiRunHead{17, 8, 2, 6, 14}
    ),
    std::make_tuple(
      String{"abcabcababc"},
      4,
      PsiRunHeadAsc{5, 1, 2, 3, 4, 0},
      PsiRunHeadInd{0, 1, 4, 5},
      PsiRunHead{11, 6, 10, 2}
    )
  )
);
