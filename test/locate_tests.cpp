//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/18/22.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/csa.h"
#include "sr-index/sr_csa.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"

using String = std::string;
using Values = std::vector<std::size_t>;

class LocateTests : public testing::TestWithParam<std::tuple<String, String, Values>> {
 protected:

  void SetUp() override {
    const auto &data = std::get<0>(GetParam());
    auto filename = sdsl::cache_file_name(key_tmp_input_, config_);
    sdsl::store_to_file(data, filename);
    register_cache_file(key_tmp_input_, config_);
  }

  void TearDown() override {
    sdsl::util::delete_all_files(config_.file_map);
  }

  sdsl::cache_config config_;
  std::string key_tmp_input_ = "data";
};

TEST_P(LocateTests, CSA) {
  sri::CSA<> index;
  sri::constructSRI<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, SRCSA) {
  sri::SrCSA<> index(6);
  sri::constructSRI<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, SRCSASlim) {
  sri::SrCSASlim<> index(6);
  sri::constructSRI<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, SRCSASlimValidMark) {
  sri::SrCSAValidMark<sri::SrCSASlim<>> index(6);
  sri::constructSRI<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, RIndex) {
  sri::RIndex<> index;
  sri::construct<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, SRIndex) {
  sri::SRIndex<> index(6);
  sri::construct<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

TEST_P(LocateTests, SRIndexValidMark) {
  sri::SRIndexValidMark<> index(6);
  sri::construct<8>(index, config_.file_map[key_tmp_input_], config_);

  const auto &pattern = std::get<1>(GetParam());
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(GetParam());
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}


INSTANTIATE_TEST_SUITE_P(
    LocateIndex,
    LocateTests,
    testing::Values(
        std::make_tuple(String{"abcabcababc\0"}, String{"ab"}, Values{6, 8, 3, 0}),
        std::make_tuple(String{"abcabcababc\0"}, String{"aba"}, Values{6}),
        std::make_tuple(String{"abcabcababc\0"}, String{"bc"}, Values{9, 4, 1})
    )
);