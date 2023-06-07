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

class BaseLocateTests : public testing::Test {
 protected:

  void Init(const String &t_data) {
    auto filename = sdsl::cache_file_name(key_tmp_input_, config_);
    sdsl::store_to_file(t_data, filename);
    register_cache_file(key_tmp_input_, config_);
  }

  void TearDown() override {
    sdsl::util::delete_all_files(config_.file_map);
  }

  sdsl::cache_config config_;
  std::string key_tmp_input_ = "data";
};

using Values = std::vector<std::size_t>;

typedef std::shared_ptr<sri::IndexBaseWithExternalStorage<>>(*TConstructor)(
    const std::string &tt_data_path, sdsl::cache_config &tt_config);

class LocateTests : public BaseLocateTests,
                    public testing::WithParamInterface<std::tuple<TConstructor, std::tuple<String, String, Values>>> {
 protected:

  void SetUp() override {
    const auto &data = std::get<0>(std::get<1>(GetParam()));
    Init(data);
    sdsl::construct_config().byte_algo_sa = sdsl::LIBDIVSUFSORT;
  }
};

TEST_P(LocateTests, Locate) {
  auto buildIndex = std::get<0>(GetParam());
  auto index = buildIndex(config_.file_map[key_tmp_input_], config_);
  const auto &info = std::get<1>(GetParam());
  const auto &pattern = std::get<1>(info);

  auto results = index->Locate(pattern);
  std::sort(results.begin(), results.end());

  auto e_results = std::get<2>(info);
  std::sort(e_results.begin(), e_results.end());
  EXPECT_EQ(results, e_results);
}

template<typename TIndex>
TConstructor createIndexBuilder() {
  return [](const std::string &tt_data_path, sdsl::cache_config &tt_config)
      -> std::shared_ptr<sri::IndexBaseWithExternalStorage<>> {
    auto index = std::make_shared<TIndex>();
    sri::construct(*index, tt_data_path, tt_config);
    return index;
  };
}

template<typename TSrIndex>
TConstructor createSrIndexBuilder() {
  return [](const std::string &tt_data_path, sdsl::cache_config &tt_config)
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
            createIndexBuilder<sri::CSA<>>(),
            createSrIndexBuilder<sri::SrCSA<>>(),
            createSrIndexBuilder<sri::SrCSAValidMark<sri::SrCSA<>>>(),
            createSrIndexBuilder<sri::SrCSAValidArea<sri::SrCSA<>>>(),
            createSrIndexBuilder<sri::SrCSASlim<>>(),
            createSrIndexBuilder<sri::SrCSAValidMark<sri::SrCSASlim<>>>(),
            createSrIndexBuilder<sri::SrCSAValidArea<sri::SrCSASlim<>>>()
        ),
        testing::Values(
            std::make_tuple(String{"abcabcababc\0"}, String{"ab"}, Values{6, 8, 3, 0}),
            std::make_tuple(String{"abcabcababc\0"}, String{"aba"}, Values{6}),
            std::make_tuple(String{"abcabcababc\0"}, String{"bc"}, Values{9, 4, 1})
        )
    )
);

template<typename TIndex>
class LocateTypedTests : public BaseLocateTests {
 public:
  void SetUp() override {
    data_ = std::make_tuple(String{"abcabcababc\0"}, String{"ab"}, Values{6, 8, 3, 0});

    const auto &data = std::get<0>(data_);
    Init(data);
    sdsl::construct_config().byte_algo_sa = sdsl::SE_SAIS;
  }

  std::tuple<String, String, Values> data_;
};

template<typename TIndex>
class RIndexLocateTypedTests : public LocateTypedTests<TIndex> {};

using RIndexes = ::testing::Types<sri::RIndex<>, sri::CSA<>>;
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
