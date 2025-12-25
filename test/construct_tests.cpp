//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 11/24/2023.
//

#include <gtest/gtest.h>

#include <sdsl/io.hpp>

#include "sr-index/config.h"
#include "sr-index/r_csa.h"
#include "sr-index/r_index.h"
#include "sr-index/sr_csa.h"

#include "base_tests.h"

using PsiRunHead = IntVector;
using PsiRunHeadAsc = IntVector;
using PsiRunHeadInd = IntVector;
using PsiRunTail = IntVector;
using PsiRunTailAsc = IntVector;
using PsiRunTailAscLink = IntVector;
using BitVector = sdsl::bit_vector;
using SDVector = sdsl::sd_vector<>;
using Marks = SDVector;
using SampleRate = std::size_t;
using SampleIdxs = SDVector;
using CumulativeRuns = IntVector;

//~~~~~~~

template <typename T>
struct Item {
  std::string key;
  T value;
  bool add_type_hash = false;
};

template <typename TIndex>
class RIndexTests : public BaseConfigTests,
                    public testing::WithParamInterface<              //
                        std::tuple<                                  //
                            typename TIndex::Alphabet::string_type,  // Input data
                            std::tuple<Item<IntVector>,              // SA
                                       Item<IntVector>,              // BWT
                                       Item<IntVector>,              // Samples
                                       Item<SDVector>                // Marks
                                       >>> {
 public:
  using Index = TIndex;
  using Data = typename TIndex::Alphabet::string_type;

 protected:
  void SetUp() override {
    const auto& data = std::get<0>(this->GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

using namespace sri::conf;

//~~~~~~~


class RIndexBytesTests : public RIndexTests<sri::RIndex<>> {
 public:
  static const sri::JSON keys;
};

const sri::JSON RIndexBytesTests::keys = sri::createDefaultKeys<Index::Alphabet::int_width>();

TEST_P(RIndexBytesTests, construct) {
  using namespace sri::conf;
  Index index;
  sri::construct(index, config_.data_path, config_);

  auto items = std::get<1>(GetParam());
  for_each_tuple(items, [this](auto&& tt_item) {
    decltype(tt_item.value) value;
    load_from_cache(value, tt_item.key, config_, tt_item.add_type_hash);

    EXPECT_THAT(value, testing::ElementsAreArray(tt_item.value)) << "Key = " << tt_item.key;
  });
}

INSTANTIATE_TEST_SUITE_P(
    Basic,
    RIndexBytesTests,
    testing::Values(std::make_tuple(  //
        RIndexBytesTests::Data{'a', 'l', 'a', 'b', 'a', 'r', 'a', 'l', 'a', 'a', 'l', 'a', 'b', 'a', 'r', 'd', 'a'},
        std::make_tuple(
            Item<IntVector> /* SA */ {
                RIndexBytesTests::keys[kSA],
                IntVector{17, 16, 8, 2, 11, 6, 0, 9, 4, 13, 3, 12, 15, 7, 1, 10, 5, 14},
            },
            Item<IntVector> /* BWT */ {
                RIndexBytesTests::keys[kBWT][kBase],
                IntVector{'a', 'd', 'l', 'l', 'l', 'r', 0, 'a', 'b', 'b', 'a', 'a', 'r', 'a', 'a', 'a', 'a', 'a'},
            },
            Item<IntVector> /* Samples */ {
                RIndexBytesTests::keys[kBWT][kTail][kTextPos],
                IntVector{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
            },
            Item<SDVector> /* Marks */ {
                RIndexBytesTests::keys[kBWT][kHead][kTextPos],
                SDVector{BitVector{0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1}},
                true,
            }))));

//~~~~~~~


class RIndexIntsTests : public RIndexTests<sri::RIndex<sri::GenericStorage, sri::Alphabet<0>>> {
 public:
  static const sri::JSON keys;
};

const sri::JSON RIndexIntsTests::keys = sri::createDefaultKeys<Index::Alphabet::int_width>();

TEST_P(RIndexIntsTests, construct) {
  using namespace sri::conf;
  Index index;
  sri::construct(index, config_.data_path, config_);

  auto items = std::get<1>(GetParam());
  for_each_tuple(items, [this](auto&& tt_item) {
    decltype(tt_item.value) value;
    load_from_cache(value, tt_item.key, config_, tt_item.add_type_hash);

    EXPECT_THAT(value, testing::ElementsAreArray(tt_item.value)) << "Key = " << tt_item.key;
  });
}

INSTANTIATE_TEST_SUITE_P(
    Basic,
    RIndexIntsTests,
    testing::Values(std::make_tuple(  //
        RIndexIntsTests::Data{'a', 'l', 'a', 'b', 'a', 'r', 'a', 'l', 'a', 'a', 'l', 'a', 'b', 'a', 'r', 'd', 'a'},
        std::make_tuple(
            Item<IntVector> /* SA */ {
                RIndexIntsTests::keys[kSA],
                IntVector{17, 16, 8, 2, 11, 6, 0, 9, 4, 13, 3, 12, 15, 7, 1, 10, 5, 14},
            },
            Item<IntVector> /* BWT */ {
                RIndexIntsTests::keys[kBWT][kBase],
                IntVector{'a', 'd', 'l', 'l', 'l', 'r', 0, 'a', 'b', 'b', 'a', 'a', 'r', 'a', 'a', 'a', 'a', 'a'},
            },
            Item<IntVector> /* Samples */ {
                RIndexBytesTests::keys[kBWT][kTail][kTextPos],
                IntVector{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
            },
            Item<SDVector> /* Marks */ {
                RIndexBytesTests::keys[kBWT][kHead][kTextPos],
                SDVector{BitVector{0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1}},
                true,
            }))));

//~~~~~~~


class BaseConstructTests : public BaseConfigTests {
 public:
  template <typename T>
  void compare(const std::string& t_key, const T& t_e_values, bool t_add_type_hash = false) const {
    T values;
    load_from_cache(values, t_key, config_, t_add_type_hash);

    EXPECT_THAT(values, testing::ElementsAreArray(t_e_values)) << "Key = " << t_key;
  }
};

class RCSATests : public BaseConstructTests,
                  public testing::WithParamInterface<
                      std::tuple<String, Psi, PsiRunHead, PsiRunTail, PsiRunTailAsc, PsiRunTailAscLink, Marks>> {
 protected:
  void SetUp() override {
    const auto& data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(RCSATests, construct) {
  using namespace sri::conf;
  sri::RCSA<> index;
  sri::constructItems(index, config_);

  compare(config_.keys[kPsi][kBase], std::get<1>(GetParam()));
  compare(config_.keys[kPsi][kHead][kTextPos], std::get<2>(GetParam()), true);
  compare(config_.keys[kPsi][kTail][kTextPos], std::get<3>(GetParam()), true);
  // compare(config_.keys[kPsi][kTail][kTextPosAsc][kIdx], std::get<4>(GetParam()));
  compare(config_.keys[kPsi][kTail][kTextPosAsc][kLink], std::get<5>(GetParam()), true);

  compare(config_.keys[kPsi][kTail][kTextPos], std::get<6>(GetParam()), true);
}

INSTANTIATE_TEST_SUITE_P(
    Basic,
    RCSATests,
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
            Marks({0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1})),
        std::make_tuple(String{"abcabcababc"},
                        Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3},
                        // Psi{4,    5,    6,    7,    8,    2,    9,    10,   11,   0,    1,    3},
                        // SA{11,    6,    8,    3,    0,    7,    9,    4,    1,    10,   5,    2},
                        // BWT{'c',  'c',  'b',  'c',  '$',  'a',  'a',  'a',  'a',  'b',  'b',  'b'},
                        PsiRunHead{11, 6, 7, 9, 10, 2},
                        PsiRunTail{11, 0, 7, 1, 5, 2},
                        PsiRunTailAsc{1, 3, 5, 4, 2, 0},
                        PsiRunTailAscLink{2, 4, 0, 5, 3, 1},
                        Marks({1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1}))));

class SRCSATests : public BaseConstructTests,
                   public testing::WithParamInterface<std::tuple<String,
                                                                 SampleRate,
                                                                 PsiRunHeadAsc,
                                                                 PsiRunHeadInd,
                                                                 PsiRunHead,
                                                                 PsiRunTail,
                                                                 PsiRunTailAscLink,
                                                                 Marks,
                                                                 SampleIdxs,
                                                                 CumulativeRuns>> {
 protected:
  void SetUp() override {
    const auto& data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(SRCSATests, construct) {
  using namespace sri::conf;

  const auto& subsample_rate = std::get<1>(GetParam());
  sri::SrCSA<> index(subsample_rate);
  sri::constructItems(index, config_);

  const auto prefix = std::to_string(subsample_rate) + "_";

  compare(config_.keys[kPsi][kHead][kTextPosAsc][kIdx], std::get<2>(GetParam()));
  compare(prefix + config_.keys[kPsi][kHead][kIdx].get<std::string>(), std::get<3>(GetParam()));
  compare(prefix + config_.keys[kPsi][kHead][kTextPos].get<std::string>(), std::get<4>(GetParam()), true);
  compare(prefix + config_.keys[kPsi][kTail][kTextPos].get<std::string>(), std::get<5>(GetParam()), true);
  compare(prefix + config_.keys[kPsi][kTail][kTextPosAsc][kLink].get<std::string>(), std::get<6>(GetParam()), true);
  compare(prefix + config_.keys[kPsi][kTail][kTextPos].get<std::string>(), std::get<7>(GetParam()), true);
  compare(prefix + config_.keys[kPsi][kHead][kIdx].get<std::string>(), std::get<8>(GetParam()), true);
  compare(config_.keys[kPsi][kCumRun].get<std::string>(), std::get<9>(GetParam()), true);
}

INSTANTIATE_TEST_SUITE_P(Basic,
                         SRCSATests,
                         testing::Values(std::make_tuple(String{"alabaralaalabarda"},
                                                         4,
                                                         PsiRunHeadAsc{3, 5, 8, 4, 7, 2, 9, 6, 1, 0},
                                                         PsiRunHeadInd{0, 1, 2, 3, 8, 9},
                                                         PsiRunHead{17, 16, 8, 2, 5, 14},
                                                         PsiRunTail{14, 17, 16, 8, 10, 5},
                                                         PsiRunTailAscLink{5, 3, 4, 0, 2, 1},
                                                         Marks({0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1}),
                                                         SampleIdxs({1, 1, 1, 1, 0, 0, 0, 0, 1, 1}),
                                                         CumulativeRuns{1, 5, 6, 7, 8, 10}),
                                         std::make_tuple(String{"abcabcababc"},
                                                         4,
                                                         PsiRunHeadAsc{5, 1, 2, 3, 4, 0},
                                                         PsiRunHeadInd{0, 1, 2, 5},
                                                         PsiRunHead{11, 6, 7, 2},
                                                         PsiRunTail{2, 11, 0, 5},
                                                         PsiRunTailAscLink{2, 0, 3, 1},
                                                         Marks({1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1}),
                                                         SampleIdxs({1, 1, 1, 0, 0, 1}),
                                                         CumulativeRuns{1, 2, 4, 6})));

using ValidMarks = sdsl::bit_vector;
using ValidAreas = IntVector;

class SRCSAValidAreaTests : public BaseConstructTests,
                            public testing::WithParamInterface<std::tuple<String, SampleRate, ValidMarks, ValidAreas>> {
 protected:
  void SetUp() override {
    const auto& data = std::get<0>(GetParam());
    Init(data, sri::SAAlgo::SDSL_LIBDIVSUFSORT);
  }
};

TEST_P(SRCSAValidAreaTests, construct) {
  using namespace sri::conf;

  const auto& subsample_rate = std::get<1>(GetParam());
  sri::SrCSAValidArea<> index(subsample_rate);
  sri::constructItems(index, config_);

  const auto prefix = std::to_string(subsample_rate) + "_";

  compare(prefix + sri::str(config_.keys[kPsi][kTail][kTextPosAsc][kValidMark]), std::get<2>(GetParam()), true);
  compare(prefix + sri::str(config_.keys[kPsi][kTail][kTextPosAsc][kValidArea]), std::get<3>(GetParam()), true);
}

INSTANTIATE_TEST_SUITE_P(
    Basic,
    SRCSAValidAreaTests,
    testing::Values(std::make_tuple(String{"alabaralaalabarda"}, 4, ValidMarks({1, 1, 1, 0, 0, 1}), ValidAreas{1, 1}),
                    std::make_tuple(String{"abcabcababc"}, 4, ValidMarks({1, 0, 1, 0}), ValidAreas{1, 4})));
