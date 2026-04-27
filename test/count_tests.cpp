//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 2026-04-26.
//

#include <fstream>

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/config.h"
#include "sr-index/r_csa.h"
#include "sr-index/r_index.h"

#include "base_tests.h"

namespace {

using DataBytes = sri::Alphabet<8>::string_type;
using DataInts = sri::Alphabet<0>::string_type;
using Range = std::pair<std::size_t, std::size_t>;

DataBytes kCorpusBytes{"abcabcababc"};
DataInts kCorpusInts{'a', 'b', 'c', 'a', 'b', 'c', 'a', 'b', 'a', 'b', 'c'};

std::vector<DataBytes> kPatternsBytes{DataBytes{"ab"}, DataBytes{"aba"}, DataBytes{"bc"}, DataBytes{"zzz"}};
std::vector<DataInts> kPatternsInts{DataInts{'a', 'b'},
                                    DataInts{'a', 'b', 'a'},
                                    DataInts{'b', 'c'},
                                    DataInts{'z', 'z', 'z'}};

std::size_t fileSize(const std::string& t_path) {
  std::ifstream in(t_path, std::ios::ate | std::ios::binary);
  return in.tellg();
}

bool prefixesEqual(const std::string& t_short, const std::string& t_long) {
  std::ifstream a(t_short, std::ios::binary);
  std::ifstream b(t_long, std::ios::binary);
  while (a) {
    char ca;
    char cb;
    if (!a.get(ca))
      break;
    if (!b.get(cb))
      return false;
    if (ca != cb)
      return false;
  }
  return true;
}

}  // namespace

class CountBytesTests : public BaseConfigTests<8> {
 protected:
  void SetUp() override {
    Init(kCorpusBytes, sri::SDSL_LIBDIVSUFSORT);
  }
};

TEST_F(CountBytesTests, RIndexCount_matches_RIndex_from_cache) {
  sri::RIndex<sri::GenericStorage, sri::Alphabet<8>> locate_index;
  sri::construct(locate_index, config_.file_map[key_tmp_input_], config_);

  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> count_index;
  // Loads the count prefix from the same cache directory.
  count_index.load(config_);

  for (const auto& pattern : kPatternsBytes) {
    Range expected = locate_index.Count(pattern);
    Range actual = count_index.Count(pattern);
    EXPECT_EQ(expected, actual) << "pattern=" << pattern;
  }
}

TEST_F(CountBytesTests, RIndexCount_loads_locate_bundle_prefix) {
  std::string locate_path = sdsl::cache_file_name("locate_bundle", config_);
  std::string count_path = sdsl::cache_file_name("count_bundle", config_);

  Range expected_ab;
  Range expected_aba;
  Range expected_bc;
  Range expected_missing;
  {
    sri::RIndex<sri::GenericStorage, sri::Alphabet<8>> locate_index;
    sri::construct(locate_index, config_.file_map[key_tmp_input_], config_);
    expected_ab = locate_index.Count(kPatternsBytes[0]);
    expected_aba = locate_index.Count(kPatternsBytes[1]);
    expected_bc = locate_index.Count(kPatternsBytes[2]);
    expected_missing = locate_index.Count(kPatternsBytes[3]);
    sdsl::store_to_file(locate_index, locate_path);
  }
  {
    sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> count_only;
    sri::construct(count_only, config_.file_map[key_tmp_input_], config_);
    sdsl::store_to_file(count_only, count_path);
  }

  // Bundle layout invariant: count bundle is byte-identical to the prefix of the locate bundle.
  EXPECT_LT(fileSize(count_path), fileSize(locate_path));
  EXPECT_TRUE(prefixesEqual(count_path, locate_path));

  // RIndexCount can read the locate bundle (consuming only the count-prefix bytes).
  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> from_locate;
  sdsl::load_from_file(from_locate, locate_path);
  EXPECT_EQ(expected_ab, from_locate.Count(kPatternsBytes[0]));
  EXPECT_EQ(expected_aba, from_locate.Count(kPatternsBytes[1]));
  EXPECT_EQ(expected_bc, from_locate.Count(kPatternsBytes[2]));
  EXPECT_EQ(expected_missing, from_locate.Count(kPatternsBytes[3]));

  // RIndexCount can read its own count-only bundle.
  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> from_count;
  sdsl::load_from_file(from_count, count_path);
  EXPECT_EQ(expected_ab, from_count.Count(kPatternsBytes[0]));
  EXPECT_EQ(expected_bc, from_count.Count(kPatternsBytes[2]));
}

TEST_F(CountBytesTests, RIndexCount_via_CountIndex_polymorphic_pointer) {
  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> count_index;
  sri::construct(count_index, config_.file_map[key_tmp_input_], config_);

  // Verifies the new CountIndex<> sibling interface dispatches Count correctly.
  sri::CountIndex<DataBytes>* iface = &count_index;
  Range r = iface->Count(kPatternsBytes[2]);
  // "bc" occurs three times in "abcabcababc". The BWT range is half-open: count = end - start.
  EXPECT_EQ(r.second - r.first, 3u);
}

//~~~~~~~


class CountIntsTests : public BaseConfigTests<0> {
 protected:
  void SetUp() override {
    Init(kCorpusInts, sri::SDSL_LIBDIVSUFSORT);
  }
};

TEST_F(CountIntsTests, RCSACount_matches_RCSA_from_cache) {
  sri::RCSA<sri::GenericStorage, sri::Alphabet<0>> locate_index;
  sri::construct(locate_index, config_.file_map[key_tmp_input_], config_);

  sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> count_index;
  count_index.load(config_);

  for (const auto& pattern : kPatternsInts) {
    Range expected = locate_index.Count(pattern);
    Range actual = count_index.Count(pattern);
    EXPECT_EQ(expected, actual);
  }
}

TEST_F(CountIntsTests, RCSACount_loads_locate_bundle_prefix) {
  std::string locate_path = sdsl::cache_file_name("locate_csa_bundle", config_);
  std::string count_path = sdsl::cache_file_name("count_csa_bundle", config_);

  std::vector<Range> expected;
  {
    sri::RCSA<sri::GenericStorage, sri::Alphabet<0>> locate_index;
    sri::construct(locate_index, config_.file_map[key_tmp_input_], config_);
    for (const auto& p : kPatternsInts)
      expected.push_back(locate_index.Count(p));
    sdsl::store_to_file(locate_index, locate_path);
  }
  {
    sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> count_only;
    sri::construct(count_only, config_.file_map[key_tmp_input_], config_);
    sdsl::store_to_file(count_only, count_path);
  }

  EXPECT_LT(fileSize(count_path), fileSize(locate_path));
  EXPECT_TRUE(prefixesEqual(count_path, locate_path));

  sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> from_locate;
  sdsl::load_from_file(from_locate, locate_path);
  for (std::size_t i = 0; i < kPatternsInts.size(); ++i) {
    EXPECT_EQ(expected[i], from_locate.Count(kPatternsInts[i]));
  }

  sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> from_count;
  sdsl::load_from_file(from_count, count_path);
  for (std::size_t i = 0; i < kPatternsInts.size(); ++i) {
    EXPECT_EQ(expected[i], from_count.Count(kPatternsInts[i]));
  }
}
