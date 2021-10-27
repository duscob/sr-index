//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/17/21.
//

#include <string>

#include <gtest/gtest.h>

#include "sr-index/toehold.h"
#include "sr-index/rle_string.hpp"
#include "sr-index/psi.h"

using BWT = std::string;
using Range = std::pair<std::size_t, std::size_t>;
using Char = char;
using Step = std::size_t;
using Data = sri::DataBackwardSearchStepBackward<Char>;
using DataFw = sri::DataBackwardSearchStepForward;

auto tie(const Data &t_data) {
  return std::tie(t_data.step, t_data.run_data.c, t_data.run_data.end);
}

auto tie(const DataFw &t_data) {
  return std::tie(t_data.step, t_data.run_data);
}

class GetLastSpecialBackwardSearchStepForPhiBackward_Test
    : public testing::TestWithParam<std::tuple<BWT, std::tuple<Range, Range, Char, Step, Data>, Data>> {
};

TEST_P(GetLastSpecialBackwardSearchStepForPhiBackward_Test, execute) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &[range, next_range, c, step, prev_data] = std::get<1>(GetParam());

  auto get_data = sri::buildComputeDataBackwardSearchStepForPhiBackward(std::cref(bwt_rle));
  auto data = get_data(range, next_range, c, step, prev_data);

  const auto &e_data = std::get<2>(GetParam());
  EXPECT_EQ(tie(data), tie(e_data));
}

INSTANTIATE_TEST_SUITE_P(
    Toehold,
    GetLastSpecialBackwardSearchStepForPhiBackward_Test,
    testing::Values(
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{0, 0}, 1, 4, Data{5, {3, 11}}), // Input
            Data{4, {1, 11}} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{1, 4}, 2, 4, Data{5, {3, 11}}), // Input
            Data{4, {2, 11}} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{5, 8}, 3, 4, Data{5, {3, 11}}), // Input
            Data{5, {3, 11}} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{0, 11}, Range{9, 11}, 4, 4, Data{5, {3, 11}}), // Input
            Data{4, {4, 11}} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{9, 11}, Range{6, 8}, 3, 1, Data{2, {4, 11}}), // Input
            Data{2, {4, 11}} // Next data
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            std::make_tuple(Range{9, 11}, Range{1, 0}, 1, 1, Data{2, {4, 11}}), // Input
            Data{2, {4, 11}} // Next data
        )
    )
);

using Psi = sdsl::int_vector<>;
using Cumulative = sdsl::int_vector<64>;
struct DataLF {
  std::size_t value = 0;

  struct Run {
    std::size_t start = 0;
    std::size_t end = 0;
  } run;

  bool operator<(const DataLF &rhs) const {
    return value < rhs.value; // || run.start < rhs.run.start || run.end < rhs.run.end || run.rank < rhs.run.rank;
  }
};

struct RunLF {
  DataLF start;
  DataLF end;
};

class GetLastSpecialBackwardSearchStepForPhiForward_Test
    : public testing::TestWithParam<
        std::tuple<std::tuple<Range, RunLF, Char, Step, DataFw>, DataFw>> {
};

TEST_P(GetLastSpecialBackwardSearchStepForPhiForward_Test, execute) {
  const auto &[range, next_range, c, step, prev_data] = std::get<0>(GetParam());

  auto get_data = sri::buildComputeDataBackwardSearchStepForPhiForward();
  auto data = get_data(range, next_range, c, step, prev_data);

  const auto &e_data = std::get<1>(GetParam());
  EXPECT_EQ(tie(data), tie(e_data));
}

INSTANTIATE_TEST_SUITE_P(
    Toehold,
    GetLastSpecialBackwardSearchStepForPhiForward_Test,
    testing::Values(
        // Symbols for Psi are in [0..sigma)
        std::make_tuple(
            std::make_tuple(Range{0, 12}, RunLF{DataLF{0, {4, 5}}, DataLF{1, {12, 12}}}, 0, 4, DataFw{5, 0}), // Input
            DataFw{4, 4} // Next data
        ),
        std::make_tuple(
            std::make_tuple(Range{0, 12}, RunLF{DataLF{1, {5, 9}}, DataLF{5, {12, 12}}}, 1, 4, DataFw{5, 0}), // Input
            DataFw{4, 5} // Next data
        ),
        std::make_tuple(
            std::make_tuple(Range{0, 12}, RunLF{DataLF{5, {2, 3}}, DataLF{9, {12, 12}}}, 2, 4, DataFw{5, 0}), // Input
            DataFw{4, 2} // Next data
        ),
        std::make_tuple(
            std::make_tuple(Range{0, 12}, RunLF{DataLF{9, {0, 2}}, DataLF{12, {12, 12}}}, 3, 4, DataFw{5, 0}), // Input
            DataFw{5, 0} // Next data
        ),
        std::make_tuple(
            std::make_tuple(Range{9, 12}, RunLF{DataLF{6, {9, 12}}, DataLF{9, {12, 12}}}, 2, 1, DataFw{2, 0}), // Input
            DataFw{2, 0} // Next data
        ),
        std::make_tuple(
            std::make_tuple(Range{9, 12}, RunLF{DataLF{0, {0, 0}}, DataLF{0, {0, 0}}}, 1, 1, DataFw{2, 0}), // Input
            DataFw{2, 0} // Next data
        )
    )
);

using SA = sdsl::int_vector<>;
using SAValue = std::size_t;

class ComputeToeholdValueForPhiBackward_Test : public testing::TestWithParam<std::tuple<BWT, SA, Data, SAValue>> {
};

TEST_P(ComputeToeholdValueForPhiBackward_Test, execute) {
  const auto &bwt = std::get<0>(GetParam());
  sri::rle_string<> bwt_rle(bwt);

  const auto &sa = std::get<1>(GetParam());
  auto get_sa = [&sa](auto tt_i) { return sa[tt_i]; };

  auto get_toehold_value = sri::buildComputeToeholdForPhiBackward(std::cref(bwt_rle), get_sa);

  const auto &data = std::get<2>(GetParam());
  auto sa_value = get_toehold_value(data);

  const auto &e_sa_value = std::get<3>(GetParam());
  EXPECT_EQ(sa_value, e_sa_value);
}

INSTANTIATE_TEST_SUITE_P(
    Toehold,
    ComputeToeholdValueForPhiBackward_Test,
    testing::Values(
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{0, {1, 11}}, // Data for pattern "1"
            11 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{0, {2, 11}}, // Data for pattern "2"
            0 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{0, {3, 11}}, // Data for pattern "3"
            1 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{0, {4, 11}}, // Data for pattern "4"
            2 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{2, {4, 11}}, // Data for pattern "2 3 4"
            0 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{0, {3, 4}}, // Data for pattern "3 2"
            7 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{1, {3, 4}}, // Data for pattern "2 3 2"
            6 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{1, {3, 11}}, // Data for pattern "2 3"
            0 // SA value
        ),
        std::make_tuple(
            BWT{4, 4, 3, 4, 1, 2, 2, 2, 2, 3, 3, 3}, // BWT
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            Data{2, {4, 11}}, // Data for pattern "2 3 4"
            0 // SA value
        )
    )
);

class ComputeToeholdValueForPhiForward_Test
    : public testing::TestWithParam<std::tuple<Psi, Cumulative, SA, DataFw, SAValue>> {
};

TEST_P(ComputeToeholdValueForPhiForward_Test, execute) {
  const auto &psi_raw = std::get<0>(GetParam());
  const auto &cum_c = std::get<1>(GetParam());
  auto psi_core = sri::PsiCoreBV(cum_c, psi_raw);

  const auto &sa = std::get<2>(GetParam());
  auto get_sa = [&sa](auto tt_i) { return sa[tt_i]; };

  auto get_toehold_value = sri::buildComputeToeholdForPhiForward(get_sa, psi_core.size());

  const auto &data = std::get<3>(GetParam());
  auto sa_value = get_toehold_value(data);

  const auto &e_sa_value = std::get<4>(GetParam());
  EXPECT_EQ(sa_value, e_sa_value);
}

INSTANTIATE_TEST_SUITE_P(
    Toehold,
    ComputeToeholdValueForPhiForward_Test,
    testing::Values(
        // Symbols for Psi are in [0..sigma)
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{0, 4}, // Data for pattern "0"
            11 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{0, 5}, // Data for pattern "1"
            6 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{0, 2}, // Data for pattern "2"
            7 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{0, 0}, // Data for pattern "3"
            10 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{2, 0}, // Data for pattern "1 2 3"
            8 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{0, 2}, // Data for pattern "2 1"
            7 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{1, 2}, // Data for pattern "1 2 1"
            6 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{1, 2}, // Data for pattern "1 2"
            6 // SA value
        ),
        std::make_tuple(
            Psi{4, 5, 6, 7, 8, 2, 9, 10, 11, 0, 1, 3}, // Psi
            Cumulative{0, 1, 5, 9, 12}, // Cumulative number of symbols
            SA{11, 6, 8, 3, 0, 7, 9, 4, 1, 10, 5, 2}, // SA
            DataFw{2, 0}, // Data for pattern "1 2 3"
            8 // SA value
        )
    )
);