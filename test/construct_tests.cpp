//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/2023.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/r_csa.h"
#include "sr-index/config.h"

#include "base_tests.h"

using PsiRunHead = IntVector;
using PsiRunTail = IntVector;
using PsiRunTailAsc = IntVector;
using PsiRunTailAscLink = IntVector;
using BitVector = sdsl::bit_vector;
using Marks = BitVector;

class ConstructRCSATests : public BaseConfigTests,
                           public testing::WithParamInterface<std::tuple<
                               String, Psi, PsiRunHead, PsiRunTail, PsiRunTailAsc, PsiRunTailAscLink, Marks
                           >> {
 protected:

  void SetUp() override {
    const auto &data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(ConstructRCSATests, construct) {
  using namespace sri::conf;
  sri::constructRCSAWithPsiRuns<8, sdsl::sd_vector<>>(config_.file_map[key_tmp_input_], config_);

  auto compare = [this](const auto &tt_key, const auto &tt_e_values) {
    sdsl::int_vector<> values;
    sdsl::load_from_cache(values, tt_key, config_);

    EXPECT_THAT(values, testing::ElementsAreArray(tt_e_values)) << "Key = " << tt_key;
  };

  auto compareBV = [this](const auto &tt_key, const auto &tt_e_values) {
    sdsl::sd_vector<> values;
    sdsl::load_from_cache(values, tt_key, config_, true);

    EXPECT_THAT(values, testing::ElementsAreArray(tt_e_values)) << "Key = " << tt_key;
  };

  compare(sdsl::conf::KEY_PSI, std::get<1>(GetParam()));
  compare(config_.keys[kPsi][kHead][kTextPos], std::get<2>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPos], std::get<3>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPosAsc][kIdx], std::get<4>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPosAsc][kLink], std::get<5>(GetParam()));

  compareBV(config_.keys[kPsi][kTail][kTextPos], std::get<6>(GetParam()));
}

INSTANTIATE_TEST_SUITE_P(
    Basic,
    ConstructRCSATests,
    testing::Values(
        std::make_tuple(
            String{"alabaralaalabarda"},
            Psi{6, 0, 7, 10, 11, 13, 14, 15, 16, 17, 8, 9, 1, 2, 3, 4, 5, 12},
//          Psi{6,    0,    7,    10,   11,   13,   14,   15,   16,   17,   8,    9,    1,    2,    3,    4,    5,    12},
//          SA{17,    16,   8,    2,    11,   6,    0,    9,    4,    13,   3,    12,   15,   7,    1,    10,   5,    14}
//          BWT{'a',  'd',  'l',  'l',  'l',  'r',  '$',  'a',  'b',  'b',  'a',  'a',  'r',  'a',  'a',  'a',  'a',  'a'},
            PsiRunHead{17, 16, 8, 2, 6, 3, 15, 7, 5, 14},
            PsiRunTail{17, 16, 8, 11, 13, 12, 15, 10, 5, 14},
            PsiRunTailAsc{8, 2, 7, 3, 5, 4, 9, 6, 1, 0},
            PsiRunTailAscLink{9, 3, 8, 4, 6, 5, 0, 7, 2, 1},
            Marks{0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}
        ),
        std::make_tuple(
            String{"abcabcababc"},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
//            Psi{4,    5,    6,    7,    8,    2,    9,    10,   11,   0,    1,    3},
//            SA{11,    6,    8,    3,    0,    7,    9,    4,    1,    10,   5,    2},
//            BWT{'c',  'c',  'b',  'c',  '$',  'a',  'a',  'a',  'a',  'b',  'b',  'b'},
            PsiRunHead{11, 6, 7, 9, 10, 2},
            PsiRunTail{11, 0, 7, 1, 5, 2},
            PsiRunTailAsc{1, 3, 5, 4, 2, 0},
            PsiRunTailAscLink{2, 4, 0, 5, 3, 1},
            Marks{1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1}
        )
    )
);