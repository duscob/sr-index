//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 1/18/22.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/r_csa.h"
#include "sr-index/sr_csa.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_index.h"
#include "sr-index/config.h"

#include "base_tests.h"

using String = std::string;

using Values = std::vector<std::size_t>;
using PatternXValues = std::tuple<String, Values>;
using ListPatternXValues = std::vector<PatternXValues>;

typedef std::shared_ptr<sri::IndexBaseWithExternalStorage<>>(*TConstructor)(
    const std::string &tt_data_path, sri::Config &tt_config);

class LocateTests : public BaseConfigTests,
                    public testing::WithParamInterface<std::tuple<
                        TConstructor,
                        std::tuple<String, ListPatternXValues>,
                        sri::SAAlgo
                    >> {
 protected:

  void SetUp() override {
    const auto &data = std::get<0>(std::get<1>(GetParam()));
    const auto &sa_algo = std::get<2>(GetParam());
#ifndef NDEBUG
    if (sa_algo == sri::SAAlgo::BIG_BWT) {
      GTEST_SKIP_("Tests with BigBWT fail in Debug mode");
    }
#endif
    Init(data, sa_algo);
  }
};

TEST_P(LocateTests, Locate) {
  auto buildIndex = std::get<0>(GetParam());
  auto index = buildIndex(config_.file_map[key_tmp_input_], config_);
  const auto &info = std::get<1>(GetParam());

  const auto &listPatternXValues = std::get<1>(info);
  for (const auto &item : listPatternXValues) {
    const auto &pattern = std::get<0>(item);

    auto results = index->Locate(pattern);
    std::sort(results.begin(), results.end());

    auto e_results = std::get<1>(item);
    std::sort(e_results.begin(), e_results.end());
    EXPECT_EQ(results, e_results) << pattern;
  }
}

template<typename TIndex>
TConstructor createIndexBuilder() {
  return [](const std::string &tt_data_path, sri::Config &tt_config)
      -> std::shared_ptr<sri::IndexBaseWithExternalStorage<>> {
    auto index = std::make_shared<TIndex>();
    sri::construct(*index, tt_data_path, tt_config);
    return index;
  };
}

template<typename TSrIndex>
TConstructor createSrIndexBuilder() {
  return [](const std::string &tt_data_path, sri::Config &tt_config)
      -> std::shared_ptr<sri::IndexBaseWithExternalStorage<>> {
    auto index = std::make_shared<TSrIndex>(6);
    sri::construct(*index, tt_data_path, tt_config);
    return index;
  };
}

INSTANTIATE_TEST_SUITE_P(
    LocateIndex,
    LocateTests,
    testing::Combine(
        testing::Values(
            createIndexBuilder<sri::RIndex<>>(),
            createSrIndexBuilder<sri::SrIndex<>>(),
            createSrIndexBuilder<sri::SrIndexValidMark<>>(),
            createSrIndexBuilder<sri::SrIndexValidArea<>>(),
            createIndexBuilder<sri::RCSA<>>(),
            createSrIndexBuilder<sri::SrCSA<>>(),
            createSrIndexBuilder<sri::SrCSAValidMark<sri::SrCSA<>>>(),
            createSrIndexBuilder<sri::SrCSAValidArea<sri::SrCSA<>>>(),
            createSrIndexBuilder<sri::SrCSASlim<>>(),
            createSrIndexBuilder<sri::SrCSAValidMark<sri::SrCSASlim<>>>(),
            createSrIndexBuilder<sri::SrCSAValidArea<sri::SrCSASlim<>>>()
        ),
        testing::Values(
            std::make_tuple(String{"abcabcababc"},
                            ListPatternXValues{
                                std::make_tuple(String{"ab"}, Values{6, 8, 3, 0}),
                                std::make_tuple(String{"aba"}, Values{6}),
                                std::make_tuple(String{"bc"}, Values{9, 4, 1})
                            }
            )
        ),
        testing::Values(
            sri::SDSL_LIBDIVSUFSORT,
            sri::BIG_BWT // Fails in Debug Mode
        )
    )
);

template<typename TIndex>
class LocateTypedTests : public BaseConfigTests {
 public:
  void SetUp() override {
    data_ = std::make_tuple(String{"abcabcababc"}, String{"ab"}, Values{6, 8, 3, 0});

    const auto &data = std::get<0>(data_);
    Init(data, sri::SDSL_SE_SAIS);
  }

  std::tuple<String, String, Values> data_;
};

template<typename TIndex>
class RIndexLocateTypedTests : public LocateTypedTests<TIndex> {};

using RIndexes = ::testing::Types<sri::RIndex<>, sri::RCSA<>>;
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

  const auto &pattern = std::get<1>(this->data_);
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(this->data_);
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

template<typename TIndex>
class SRIndexLocateTypedTests : public LocateTypedTests<TIndex> {};

using SRIndexes = ::testing::Types<sri::SrIndex<>,
                                   sri::SrIndexValidMark<>,
                                   sri::SrIndexValidArea<>,
                                   sri::SrCSA<>,
                                   sri::SrCSAValidMark<sri::SrCSA<>>,
                                   sri::SrCSAValidArea<sri::SrCSA<>>,
                                   sri::SrCSASlim<>,
                                   sri::SrCSAValidMark<sri::SrCSASlim<>>,
                                   sri::SrCSAValidArea<sri::SrCSASlim<>>>;
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

  const auto &pattern = std::get<1>(this->data_);
  auto results = index.Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(this->data_);
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}
