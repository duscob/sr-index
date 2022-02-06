//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 10/21/21.
//

#ifndef SRI_TEST_PSI_BASE_TESTS_H_
#define SRI_TEST_PSI_BASE_TESTS_H_

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/util.hpp>
#include <sdsl/config.hpp>
#include <sdsl/io.hpp>

using BWT = sdsl::int_vector<8>;

class BaseAlphabetTests : public testing::Test {
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
using Range = std::pair<std::size_t, std::size_t>;
using Char = unsigned char;

#endif //SRI_TEST_PSI_BASE_TESTS_H_
