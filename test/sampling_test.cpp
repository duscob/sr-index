//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/21.
//

#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <sdsl/int_vector.hpp>

#include "sr-index/sampling.h"

using SamplingSize = std::size_t;
using Values = std::vector<uint32_t>;
//using Indexes = std::vector<uint32_t>;
using Indexes = sdsl::int_vector<>;

class ComputeSampling_Test
    : public testing::TestWithParam<std::tuple<SamplingSize, Values, Indexes, Indexes, Indexes>> {
};

class ComputeSamplingForward_Test : public ComputeSampling_Test {};

TEST_P(ComputeSamplingForward_Test, execute) {
  const auto &s = std::get<0>(GetParam());
  const auto &values = std::get<1>(GetParam());
  const auto &idxs = std::get<2>(GetParam());
  const auto &req_idxs = std::get<3>(GetParam());

  auto sampled_idxs = sri::computeSampling(s, idxs.begin(), idxs.end(), values, req_idxs);

  const auto &e_sampled_idxs = std::get<4>(GetParam());
  EXPECT_THAT(sampled_idxs, testing::ElementsAreArray(e_sampled_idxs));
}

INSTANTIATE_TEST_SUITE_P(
    Sampling,
    ComputeSamplingForward_Test,
    testing::Values(
        // Text: abcabcababc$   ---   Samples: BWT-run ends
        std::make_tuple(1, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{}, Indexes{4, 5, 2, 0, 1, 3}),
        std::make_tuple(2, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{}, Indexes{4, 2, 0, 1, 3}),
        std::make_tuple(4, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{}, Indexes{4, 2, 0, 1, 3}),
        std::make_tuple(6, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{}, Indexes{4, 0, 3}),
        std::make_tuple(8, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{}, Indexes{4, 1, 3}),
        std::make_tuple(6, Values{5, 7, 2, 11, 0, 1}, Indexes{4, 5, 2, 0, 1, 3}, Indexes{1, 2}, Indexes{4, 2, 1, 3}),

        // Text: alabaralaalabarda$   ---   Samples: BWT-run ends
        std::make_tuple(1,
                        Values{16, 15, 10, 5, 17, 8, 12, 11, 14, 13}, // Samples (BWT-runs ends)
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{},
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4}), // Expected sub-sampled indexes
        std::make_tuple(4,
                        Values{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4},
                        Indexes{},
                        Indexes{3, 5, 6, 0, 4}),
        std::make_tuple(8,
                        Values{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4},
                        Indexes{},
                        Indexes{3, 9, 4}),
        std::make_tuple(4,
                        Values{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4},
                        Indexes{6, 3},
                        Indexes{3, 5, 6, 0, 4}),
        std::make_tuple(8,
                        Values{16, 15, 10, 5, 17, 8, 12, 11, 14, 13},
                        Indexes{3, 5, 2, 7, 6, 9, 8, 1, 0, 4},
                        Indexes{6, 3},
                        Indexes{3, 6, 4})
    )
);

class ComputeSamplingBackward_Test : public ComputeSampling_Test {};

TEST_P(ComputeSamplingBackward_Test, execute) {
  const auto &s = std::get<0>(GetParam());
  const auto &values = std::get<1>(GetParam());
  auto idxs = std::get<2>(GetParam());
  const auto &req_idxs = std::get<3>(GetParam());

//  auto sampled_idxs = sri::computeSampling(s, idxs.rbegin(), idxs.rend(), values, req_idxs);
  auto sampled_idxs = sri::computeSampling(
      s, std::reverse_iterator(idxs.end()), std::reverse_iterator(idxs.begin()), values, req_idxs);

  const auto &e_sampled_idxs = std::get<4>(GetParam());
  EXPECT_THAT(sampled_idxs, testing::ElementsAreArray(e_sampled_idxs));
}

INSTANTIATE_TEST_SUITE_P(
    Sampling,
    ComputeSamplingBackward_Test,
    testing::Values(
        // Text: abcabcababc$    ---   Samples: BWT-run starts
        std::make_tuple(1, Values{10, 7, 2, 11, 6, 9}, Indexes{2, 4, 1, 5, 0, 3}, Indexes{}, Indexes{3, 0, 5, 1, 4, 2}),
        std::make_tuple(2, Values{10, 7, 2, 11, 6, 9}, Indexes{2, 4, 1, 5, 0, 3}, Indexes{}, Indexes{3, 5, 1, 4, 2}),
        std::make_tuple(4, Values{10, 7, 2, 11, 6, 9}, Indexes{2, 4, 1, 5, 0, 3}, Indexes{}, Indexes{3, 1, 4, 2}),
        std::make_tuple(8, Values{10, 7, 2, 11, 6, 9}, Indexes{2, 4, 1, 5, 0, 3}, Indexes{}, Indexes{3, 4, 2}),
        std::make_tuple(8, Values{10, 7, 2, 11, 6, 9}, Indexes{2, 4, 1, 5, 0, 3}, Indexes{1, 4}, Indexes{3, 1, 4, 2}),

        // Text: alabaralaalabarda$   ---   Samples: BWT-run starts
        std::make_tuple(1,
                        Values{16, 15, 7, 5, 17, 8, 3, 2, 14, 6}, // Samples (BWT-run starts)
                        Indexes{7, 6, 3, 9, 2, 5, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{},
                        Indexes{4, 0, 1, 8, 5, 2, 9, 3, 6, 7}), // Expected sub-sampled indexes
        std::make_tuple(4,
                        Values{16, 15, 7, 5, 17, 8, 3, 2, 14, 6}, // Samples (BWT-run starts)
                        Indexes{7, 6, 3, 9, 2, 5, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{},
                        Indexes{4, 8, 5, 3, 7}), // Expected sub-sampled indexes
        std::make_tuple(8,
                        Values{16, 15, 7, 5, 17, 8, 3, 2, 14, 6}, // Samples (BWT-run starts)
                        Indexes{7, 6, 3, 9, 2, 5, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{},
                        Indexes{4, 8, 9, 7}), // Expected sub-sampled indexes
        std::make_tuple(4,
                        Values{16, 15, 7, 5, 17, 8, 3, 2, 14, 6}, // Samples (BWT-run starts)
                        Indexes{7, 6, 3, 9, 2, 5, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{8, 0},
                        Indexes{4, 0, 8, 5, 3, 7}), // Expected sub-sampled indexes
        std::make_tuple(8,
                        Values{16, 15, 7, 5, 17, 8, 3, 2, 14, 6}, // Samples (BWT-run starts)
                        Indexes{7, 6, 3, 9, 2, 5, 8, 1, 0, 4}, // Sample indexes sorted
                        Indexes{8, 0},
                        Indexes{4, 0, 8, 9, 7}) // Expected sub-sampled indexes
    )
);