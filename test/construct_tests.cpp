//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/2023.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/r_csa.h"
#include "sr-index/sr_csa.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"
#include "sr-index/config.h"

#include "base_tests.h"

class ConstructRCSATests : public BaseConfigTests,
                           public testing::WithParamInterface<std::tuple<
                               String,
                               Psi
                           >> {
 protected:

  void SetUp() override {
    const auto &data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(ConstructRCSATests, construct) {
  sri::constructRCSA<8, sdsl::sd_vector<>>(config_.file_map[key_tmp_input_], config_);

  auto compare = [this](const auto &tt_key, const auto &tt_e_values) {
    sdsl::int_vector<> values;
    sdsl::load_from_cache(values, tt_key, config_);

    EXPECT_THAT(values, testing::ElementsAreArray(tt_e_values)) << tt_key;
  };

  compare(sdsl::conf::KEY_PSI, std::get<1>(GetParam()));
}

INSTANTIATE_TEST_SUITE_P(
    LocateIndex,
    ConstructRCSATests,
    testing::Values(
        std::make_tuple(
            String{"alabaralaalabarda"},
            Psi{6, 0, 7, 10, 11, 13, 14, 15, 16, 17, 8, 9, 1, 2, 3, 4, 5, 12}
//            BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
        ),
        std::make_tuple(
            String{"abcabcababc"},
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}
//            BWT{'c', 'c', 'b', 'c', '$', 'a', 'a', 'a', 'a', 'b', 'b', 'b'},
        )
    )
);