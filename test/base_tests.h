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

#include "sr-index/config.h"

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

using IntVector = sdsl::int_vector<>;
using Psi = IntVector;
using Range = std::pair<std::size_t, std::size_t>;
using Char = unsigned char;
using String = std::string;

class BaseConfigTests : public testing::Test {
 protected:

  void Init(const String &t_data, sri::SAAlgo t_sa_algo) {
    auto filename = sdsl::cache_file_name(key_tmp_input_, config_);
    sdsl::store_to_file(t_data, filename);
    register_cache_file(key_tmp_input_, config_);

    config_.data_path = filename;
    config_.sa_algo = t_sa_algo;
  }

  void TearDown() override {
    sdsl::util::delete_all_files(config_.file_map);
  }

  sri::Config config_;
  std::string key_tmp_input_ = "data";
};

#endif //SRI_TEST_PSI_BASE_TESTS_H_
