//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/18/22.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/config.h"
#include "sr-index/r_csa.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_csa.h"
#include "sr-index/sr_csa_bwt.h"
#include "sr-index/sr_idx_generic.h"
#include "sr-index/sr_index.h"

#include "base_tests.h"

using String = std::string;

using Values = std::vector<std::size_t>;
template <typename TData>
using PatternXValues = std::tuple<TData, Values>;
template <typename TData>
using ListPatternXValues = std::vector<PatternXValues<TData>>;

template <typename TData>
using TConstructor = std::function<std::shared_ptr<sri::LocateIndex<TData>>(const std::string&, sri::Config&)>;

template <typename TAlphabet>
class LocateTests
    : public BaseConfigTests<TAlphabet::int_width>,
      public testing::WithParamInterface<
          std::tuple<TConstructor<typename TAlphabet::string_type>,
                     std::tuple<typename TAlphabet::string_type, ListPatternXValues<typename TAlphabet::string_type>>,
                     sri::SAAlgo>> {
 protected:
  void SetUp() override {
    const auto& data = std::get<0>(std::get<1>(this->GetParam()));
    const auto& sa_algo = std::get<2>(this->GetParam());
#ifndef NDEBUG
    if (sa_algo == sri::SAAlgo::BIG_BWT) {
      GTEST_SKIP_("Tests with BigBWT fail in Debug mode");
    }
#endif
    this->Init(data, sa_algo);
  }
};

//~~~~~~~


template <typename TData, typename TIndex>
TConstructor<TData> createIndexBuilder() {
  return [](const std::string& tt_data_path, sri::Config& tt_config) -> std::shared_ptr<sri::LocateIndex<TData>> {
    auto index = std::make_shared<TIndex>();
    sri::construct(*index, tt_data_path, tt_config);
    return index;
  };
}

template <typename TData, typename TSrIndex>
TConstructor<TData> createSrIndexBuilder() {
  return [](const std::string& tt_data_path, sri::Config& tt_config) -> std::shared_ptr<sri::LocateIndex<TData>> {
    auto index = std::make_shared<TSrIndex>(6);
    sri::construct(*index, tt_data_path, tt_config);
    return index;
  };
}

//~~~~~~~


using DataBytes = sri::Alphabet<8>::string_type;

class LocateBytesTests : public LocateTests<sri::Alphabet<8>> {};

TEST_P(LocateBytesTests, Locate) {
  auto buildIndex = std::get<0>(GetParam());
  auto index = buildIndex(config_.file_map[key_tmp_input_], config_);
  const auto& info = std::get<1>(GetParam());

  const auto& listPatternXValues = std::get<1>(info);
  for (const auto& item : listPatternXValues) {
    const auto& pattern = std::get<0>(item);

    auto results = index->Locate(pattern);
    std::sort(results.begin(), results.end());

    auto e_results = std::get<1>(item);
    std::sort(e_results.begin(), e_results.end());
    EXPECT_EQ(results, e_results) << pattern;
  }
}

INSTANTIATE_TEST_SUITE_P(
    LocateIndex,
    LocateBytesTests,
    testing::Combine(
        testing::Values(createIndexBuilder<DataBytes, sri::RIndex<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrIndex<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrIndexValidMark<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrIndexValidArea<>>(),
                        createIndexBuilder<DataBytes, sri::RCSABWTRun<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRun<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRunValidMark<sri::SrCSABWTRun<>>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRunValidArea<sri::SrCSABWTRun<>>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRunSlim<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRunValidMark<sri::SrCSABWTRunSlim<>>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSABWTRunValidArea<sri::SrCSABWTRunSlim<>>>(),
                        createIndexBuilder<DataBytes, sri::RCSA<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSA<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSAValidMark<>>(),
                        createSrIndexBuilder<DataBytes, sri::SrCSAValidArea<>>()),
        testing::Values(std::make_tuple(DataBytes{"abcabcababc"},
                                        ListPatternXValues<DataBytes>{
                                            std::make_tuple(DataBytes{"ab"}, Values{6, 8, 3, 0}),
                                            std::make_tuple(DataBytes{"aba"}, Values{6}),
                                            std::make_tuple(DataBytes{"bc"}, Values{9, 4, 1}),
                                        })),
        testing::Values(sri::SDSL_LIBDIVSUFSORT,
                        sri::BIG_BWT  // Fails in Debug Mode
                        )));

//~~~~~~~


using DataInts = sri::Alphabet<0>::string_type;

class LocateIntsTests : public LocateTests<sri::Alphabet<0>> {};

TEST_P(LocateIntsTests, Locate) {
  auto buildIndex = std::get<0>(GetParam());
  auto index = buildIndex(config_.file_map[key_tmp_input_], config_);
  const auto& info = std::get<1>(GetParam());

  const auto& listPatternXValues = std::get<1>(info);
  for (const auto& item : listPatternXValues) {
    const auto& pattern = std::get<0>(item);

    auto results = index->Locate(pattern);
    std::sort(results.begin(), results.end());

    auto e_results = std::get<1>(item);
    std::sort(e_results.begin(), e_results.end());
    EXPECT_EQ(results, e_results);
  }
}

INSTANTIATE_TEST_SUITE_P(
    LocateIndex,
    LocateIntsTests,
    testing::Combine(
        testing::Values(createIndexBuilder<DataInts, sri::RIndex<sri::GenericStorage, sri::Alphabet<0>>>(),
                        createSrIndexBuilder<DataInts, sri::SrIndex<sri::GenericStorage, sri::Alphabet<0>>>(),
                        createSrIndexBuilder<DataInts, sri::SrIndexValidMark<sri::GenericStorage, sri::Alphabet<0>>>(),
                        createSrIndexBuilder<DataInts, sri::SrIndexValidArea<sri::GenericStorage, sri::Alphabet<0>>>()
                        // createIndexBuilder<String, sri::RCSABWTRun<>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRun<>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRunValidMark<sri::SrCSABWTRun<>>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRunValidArea<sri::SrCSABWTRun<>>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRunSlim<>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRunValidMark<sri::SrCSABWTRunSlim<>>>(),
                        // createSrIndexBuilder<String, sri::SrCSABWTRunValidArea<sri::SrCSABWTRunSlim<>>>(),
                        // createIndexBuilder<DataInts, sri::RCSA<sri::GenericStorage, sri::Alphabet<0>>>()
                        // createSrIndexBuilder<String, sri::SrCSA<>>(),
                        // createSrIndexBuilder<String, sri::SrCSAValidMark<>>(),
                        // createSrIndexBuilder<IntVector, sri::SrCSAValidArea<>>()
                        ),
        testing::Values(std::make_tuple(DataInts{'a', 'b', 'c', 'a', 'b', 'c', 'a', 'b', 'a', 'b', 'c'},
                                        ListPatternXValues<DataInts>{
                                            std::make_tuple(DataInts{'a', 'b'}, Values{6, 8, 3, 0}),
                                            std::make_tuple(DataInts{'a', 'b', 'a'}, Values{6}),
                                            std::make_tuple(DataInts{'b', 'c'}, Values{9, 4, 1}),
                                        })),
        testing::Values(sri::SDSL_LIBDIVSUFSORT)));

//~~~~~~~


template <typename TIndex>
class LocateTypedTests : public BaseConfigTests<8> {
 public:
  void SetUp() override {
    data_ = std::make_tuple(String{"abcabcababc"}, String{"ab"}, Values{6, 8, 3, 0});

    const auto& data = std::get<0>(data_);
    Init(data, sri::SDSL_SE_SAIS);
  }

  std::tuple<String, String, Values> data_;
};

template <typename TIndex>
class RIndexLocateTypedTests : public LocateTypedTests<TIndex> {};

using RIndexes = ::testing::Types<sri::RIndex<>,
                                  sri::RCSABWTRun<>,
                                  sri::RCSA<>,
                                  sri::SrIdxGeneric<sri::SrIndex<>, 2>,
                                  sri::SrIdxGeneric<sri::SrCSABWTRun<>, 4>>;
TYPED_TEST_SUITE(RIndexLocateTypedTests, RIndexes);

TYPED_TEST(RIndexLocateTypedTests, serialize) {
  auto key_index = "index";
  {
    TypeParam index;
    sri::construct(index, this->config_.file_map[this->key_tmp_input_], this->config_);
    sdsl::store_to_cache(index, key_index, this->config_);
  }

  TypeParam index;
  sdsl::load_from_cache(index, key_index, this->config_);

  const auto& pattern = std::get<1>(this->data_);
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(this->data_);
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

//~~~~~~~


template <typename TIndex>
class SRIndexLocateTypedTests : public LocateTypedTests<TIndex> {};

using SRIndexes = ::testing::Types<sri::SrIndex<>,
                                   sri::SrIndexValidMark<>,
                                   sri::SrIndexValidArea<>,
                                   sri::SrCSABWTRun<>,
                                   sri::SrCSABWTRunValidMark<sri::SrCSABWTRun<>>,
                                   sri::SrCSABWTRunValidArea<sri::SrCSABWTRun<>>,
                                   sri::SrCSABWTRunSlim<>,
                                   sri::SrCSABWTRunValidMark<sri::SrCSABWTRunSlim<>>,
                                   sri::SrCSABWTRunValidArea<sri::SrCSABWTRunSlim<>>,
                                   sri::SrCSA<>,
                                   sri::SrCSAValidMark<>,
                                   sri::SrCSAValidArea<>>;
TYPED_TEST_SUITE(SRIndexLocateTypedTests, SRIndexes);

TYPED_TEST(SRIndexLocateTypedTests, serialize) {
  auto key_index = "index";
  {
    TypeParam index(6);
    sri::construct(index, this->config_.file_map[this->key_tmp_input_], this->config_);
    sdsl::store_to_cache(index, key_index, this->config_);
  }

  TypeParam index;
  sdsl::load_from_cache(index, key_index, this->config_);

  const auto& pattern = std::get<1>(this->data_);
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(this->data_);
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}
