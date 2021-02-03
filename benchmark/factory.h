//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/28/20.
//

#ifndef RI_BENCHMARK_FACTORY_H_
#define RI_BENCHMARK_FACTORY_H_

#include <iostream>

#include <sdsl/config.hpp>

#include <ri/rle_string.hpp>
#include <ri/r_index.h>
#include <ri/tools.h>
#include <ri/bwt.h>
#include <ri/predecessor.h>
#include <ri/phi.h>

#include "definitions.h"

class Factory {
 public:

  enum class IndexEnum {
    RIndex = 0,
    RIndexSampled,
    RIndexSampledWithTrustedMarks,
    RIndexSampledWithTrustedAreas
  };

  Factory(const sdsl::cache_config &t_config) : config_{t_config} {
    std::ofstream breakdown_csv("breakdown.csv");
    breakdown_csv << "BWT" << std::endl;

    //Loading run-length encoded BWT
    load(bwt_rle_, ri::KEY_BWT_RLE);
    seq_size_ = bwt_rle_.item.size();

    // Loading F array
    load(f_, ri::KEY_F);

    // Loading original r-index components
    {
      auto &components = r_index_packs_[0];
      load(components.tails_in_text, ri::KEY_BWT_TAILS_TEXT_POS);

      load(components.heads_in_text_bv, ri::KEY_BWT_HEADS_TEXT_POS + "_bv_sd");
//      load(components.heads_in_text_bv, ri::KEY_BWT_HEADS_TEXT_POS + "_bv");

      components.rank_heads_in_text_bv.item =
          decltype(components.rank_heads_in_text_bv.item)(&components.heads_in_text_bv.item);
      components.rank_heads_in_text_bv.computeSize();
      components.select_heads_in_text_bv.item =
          decltype(components.select_heads_in_text_bv.item)(&components.heads_in_text_bv.item);
      components.select_heads_in_text_bv.computeSize();

      load(components.tail_idxs_by_heads_in_text, ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT);
    }

    // Loading sampled r-index components
    for (int i = 4; i <= 512; i *= 2) {
      auto &components = r_index_packs_[i];

      auto prefix = std::to_string(i) + "_";

      load(components.tails_in_text, prefix + ri::KEY_BWT_TAILS_TEXT_POS_SAMPLED);

      load(components.heads_in_text_bv, prefix + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv_sd");
//      load(components.heads_in_text_bv, prefix + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv");
      components.rank_heads_in_text_bv.item =
          decltype(components.rank_heads_in_text_bv.item)(&components.heads_in_text_bv.item);
      components.rank_heads_in_text_bv.computeSize();
      components.select_heads_in_text_bv.item =
          decltype(components.select_heads_in_text_bv.item)(&components.heads_in_text_bv.item);
      components.select_heads_in_text_bv.computeSize();

      load(components.tail_idxs_by_heads_in_text, prefix + ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT);

      load(components.sampled_tails_idx_bv, prefix + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv_sd");
//      load(components.sampled_tails_idx_bv, prefix + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv");
      components.rank_sampled_tails_idx_bv.item =
          decltype(components.rank_sampled_tails_idx_bv.item)(&components.sampled_tails_idx_bv.item);
      components.rank_sampled_tails_idx_bv.computeSize();

//      load(components.marked_sampled_idxs_bv, prefix + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv_sd");
      load(components.marked_sampled_idxs_bv, prefix + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv");
      components.rank_marked_sampled_idxs_bv.item =
          decltype(components.rank_marked_sampled_idxs_bv.item)(&components.marked_sampled_idxs_bv.item);

      load(components.head_marked_sample_trusted_areas, prefix + ri::KEY_BWT_HEADS_MARKED_SAMPLED_TRUSTED_AREA_IN_TEXT);

    }
  }

  auto SequenceSize() const {
    return seq_size_;
  }

  struct Config {
    IndexEnum index;
    std::size_t sampling_size;
  };

  std::pair<std::shared_ptr<ri::LocateIndex>, std::size_t> make(const Config &t_config) const {
    return internal_make(t_config);
  }

 private:
  struct RIndexComponents;

  auto makeLF() const {
    auto get_rank_of_char = ri::buildRankOfChar(std::cref(bwt_rle_.item));
    auto get_f = ri::buildRandomAccessForContainer(std::cref(f_.item));
    auto max_c = f_.item.size() - 1;

    return ri::buildLF(get_rank_of_char, get_f, seq_size_, max_c);
  }

  auto makeGetSampleForSAPosition(std::size_t t_s) const {
    const auto &components = r_index_packs_.at(t_s);

    auto get_run_of_sa_pos = ri::buildRunOfSAPosition(std::cref(bwt_rle_.item));

    auto get_is_run_sampled = ri::buildRandomAccessForContainer(std::cref(components.sampled_tails_idx_bv.item));

    auto get_sample = ri::buildRandomAccessForContainer(std::cref(components.tails_in_text.item));
    auto get_sample_for_bwt_run =
        ri::buildGetSampleForBWTRun(std::cref(components.rank_sampled_tails_idx_bv.item), get_sample);

    return ri::buildGetSampleForSAPosition(get_run_of_sa_pos, get_is_run_sampled, get_sample_for_bwt_run);

  }

  auto makeGetLastValue() const {
    auto get_run_of_sa_pos = ri::buildRunOfSAPosition(std::cref(bwt_rle_.item));
    auto always_is_sampled = [](const auto &p) { return true; };
    auto get_sample_for_bwt = ri::buildRandomAccessForContainer(std::cref(r_index_packs_.at(0).tails_in_text.item));
    auto get_sample_for_sa_pos =
        ri::buildGetSampleForSAPosition(get_run_of_sa_pos, always_is_sampled, get_sample_for_bwt);

    return ri::buildGetLastValue(std::cref(bwt_rle_.item), get_sample_for_sa_pos);
  }

  auto makeGetLastValue(std::size_t t_s) const {
    return ri::buildGetLastValue(std::cref(bwt_rle_.item), makeGetSampleForSAPosition(t_s));
  }

  template<typename TGetPredToRun, typename TSampledTailValidator>
  auto makePhi(const RIndexComponents &t_components,
               const TGetPredToRun &t_get_pred_to_run,
               const TSampledTailValidator &t_sampled_tail_validator) const {
    auto predecessor = ri::buildCircularPredecessor(std::cref(t_components.rank_heads_in_text_bv.item),
                                                    std::cref(t_components.select_heads_in_text_bv.item),
                                                    t_components.heads_in_text_bv.item.size());
    auto get_sample = ri::buildRandomAccessForContainer(std::cref(t_components.tails_in_text.item));

    return ri::buildPhi(predecessor, t_get_pred_to_run, get_sample, t_sampled_tail_validator, seq_size_);
  }

  template<typename TGetPredToRun>
  auto makePhi(const RIndexComponents &t_components, const TGetPredToRun &t_get_pred_to_run) const {
    ri::SampleValidatorDefault sampled_tail_validator_default;
    return makePhi(t_components, t_get_pred_to_run, sampled_tail_validator_default);
  }

  auto makePhi() const {
    const auto &components = r_index_packs_.at(0);

    auto get_pred_to_run = ri::buildRandomAccessForTwoContainersDefault(
        std::cref(components.tail_idxs_by_heads_in_text.item), true);

    return makePhi(components, get_pred_to_run);
  }

  template<typename TPhi>
  auto makePhiForRange(std::size_t t_s, const TPhi &t_phi) const {
    // Split in runs
    auto split_in_runs = ri::buildSplitInRuns(std::cref(bwt_rle_.item));

    return ri::buildPhiForRange(t_phi, split_in_runs, makeLF(), makeGetSampleForSAPosition(t_s), t_s, seq_size_);
  }

  auto makeComputeAllValuesWithPhi() const {
    return ri::buildComputeAllValuesWithPhi(makePhi());
  }

  template<typename TPhiForRange>
  auto makeComputeAllValuesWithPhiForRange(std::size_t t_s, const TPhiForRange &t_phi_for_range) const {
    auto get_char = ri::buildRandomAccessForContainer(std::cref(bwt_rle_.item));
    auto get_rank_of_char = ri::buildRankOfChar(std::cref(bwt_rle_.item));
    auto get_f = ri::buildRandomAccessForContainer(std::cref(f_.item));
    auto lf = ri::buildBasicLF(get_char, get_rank_of_char, get_f);

    auto get_value_for_sa_pos = ri::buildGetValueForSAPosition(makeGetSampleForSAPosition(t_s), lf, seq_size_);
    auto compute_final_value = ri::buildComputeFinalValueWithLastSampledValue(get_value_for_sa_pos);

    return ri::buildComputeAllValuesWithPhiForRange(t_phi_for_range, compute_final_value);
  }

  auto sizeBasicComponents() const {
    return bwt_rle_.size_in_bytes + f_.size_in_bytes;
  }

  auto sizeRIndexComponents(const RIndexComponents &t_components) const {
    return t_components.tails_in_text.size_in_bytes
        + t_components.heads_in_text_bv.size_in_bytes
        + t_components.rank_heads_in_text_bv.size_in_bytes
        + t_components.select_heads_in_text_bv.size_in_bytes
        + t_components.tail_idxs_by_heads_in_text.size_in_bytes
        + t_components.sampled_tails_idx_bv.size_in_bytes
        + t_components.rank_sampled_tails_idx_bv.size_in_bytes;
  }

  auto sizeRIndexComponentsWithTrustedMarks(const RIndexComponents &t_components) const {
    return sizeRIndexComponents(t_components) + t_components.marked_sampled_idxs_bv.size_in_bytes;
  }

  auto sizeRIndexComponentsWithTrustedAreas(const RIndexComponents &t_components) const {
    return sizeRIndexComponents(t_components)
        + t_components.marked_sampled_idxs_bv.size_in_bytes
        + t_components.rank_marked_sampled_idxs_bv.size_in_bytes
        + t_components.head_marked_sample_trusted_areas.size_in_bytes;
  }

  std::pair<std::shared_ptr<ri::LocateIndex>, std::size_t> internal_make(const Config &t_config) const {
    switch (t_config.index) {
      case IndexEnum::RIndex: {
        const auto &components = r_index_packs_.at(0);

        auto sa_end_value =
            ri::GetOptionalValue(components.tails_in_text.item[components.tails_in_text.item.size() - 1] + 1);

        return {ri::buildSharedPtrRIndex(makeLF(),
                                         makeGetLastValue(),
                                         makeComputeAllValuesWithPhi(),
                                         seq_size_,
                                         sa_end_value),
                sizeBasicComponents() + sizeRIndexComponents(components)};
      }

      case IndexEnum::RIndexSampled: {
        auto s = t_config.sampling_size;
        const auto &components = r_index_packs_.at(s);

        auto get_pred_to_run = ri::buildRandomAccessForTwoContainersDefault(
            std::cref(components.tail_idxs_by_heads_in_text.item), false);
        auto phi = makePhi(components, get_pred_to_run);

        auto sa_end_value =
            ri::GetOptionalValue(components.tails_in_text.item[components.tails_in_text.item.size() - 1] + 1);

        return {ri::buildSharedPtrRIndex(makeLF(),
                                         makeGetLastValue(s),
                                         makeComputeAllValuesWithPhiForRange(s, makePhiForRange(s, phi)),
                                         seq_size_,
                                         sa_end_value),
                sizeBasicComponents() + sizeRIndexComponents(components)};
      }

      case IndexEnum::RIndexSampledWithTrustedMarks: {
        auto s = t_config.sampling_size;
        const auto &components = r_index_packs_.at(s);

        auto get_pred_to_run = ri::buildRandomAccessForTwoContainers(
            std::cref(components.tail_idxs_by_heads_in_text.item), std::cref(components.marked_sampled_idxs_bv.item));
        auto phi = makePhi(components, get_pred_to_run);

        auto sa_end_value =
            ri::GetOptionalValue(components.tails_in_text.item[components.tails_in_text.item.size() - 1] + 1);

        return {ri::buildSharedPtrRIndex(makeLF(),
                                         makeGetLastValue(s),
                                         makeComputeAllValuesWithPhiForRange(s, makePhiForRange(s, phi)),
                                         seq_size_,
                                         sa_end_value),
                sizeBasicComponents() + sizeRIndexComponentsWithTrustedMarks(components)};
      }

      case IndexEnum::RIndexSampledWithTrustedAreas: {
        auto s = t_config.sampling_size;
        const auto &components = r_index_packs_.at(s);

        auto get_pred_to_run = ri::buildRandomAccessForTwoContainers(
            std::cref(components.tail_idxs_by_heads_in_text.item), std::cref(components.marked_sampled_idxs_bv.item));
        auto sample_validator = ri::buildSampleValidator(
            std::cref(components.rank_marked_sampled_idxs_bv.item),
            ri::buildRandomAccessForContainer(std::cref(components.head_marked_sample_trusted_areas.item)));
        auto phi = makePhi(components, get_pred_to_run, sample_validator);

        auto sa_end_value =
            ri::GetOptionalValue(components.tails_in_text.item[components.tails_in_text.item.size() - 1] + 1);

        return {ri::buildSharedPtrRIndex(makeLF(),
                                         makeGetLastValue(s),
                                         makeComputeAllValuesWithPhiForRange(s, makePhiForRange(s, phi)),
                                         seq_size_,
                                         sa_end_value),
                sizeBasicComponents() + sizeRIndexComponentsWithTrustedAreas(components)};
      }
    }

    exit(4);
  }

  //**********
  //**********
  //**********
  sdsl::cache_config config_;

  template<typename T>
  struct Item {
    T item;
    std::size_t size_in_bytes = 0;

    void computeSize() {
      size_in_bytes = sdsl::size_in_bytes(item);
    }
  };

  template<typename T>
  void load(Item<T> &t_item, const std::string &t_key) {
    if (!sdsl::cache_file_exists(t_key, config_))
      std::cerr << "ERROR: File '" << sdsl::cache_file_name(t_key, config_) << "' not exist!!!";
    sdsl::load_from_cache(t_item.item, t_key, config_);
    t_item.size_in_bytes = sdsl::size_in_bytes(t_item.item);
  }

  std::size_t seq_size_;

  Item<ri::rle_string<>> bwt_rle_; // Run-length encoded BWT

  Item<std::vector<std::size_t>> f_; // F Array

  struct RIndexComponents {
    Item<sdsl::int_vector<>> tails_in_text;

    using BV1 = sdsl::sd_vector<>;
//    using BV1 = sdsl::bit_vector;
    Item<BV1> heads_in_text_bv;
    Item<BV1::rank_1_type> rank_heads_in_text_bv;
    Item<BV1::select_1_type> select_heads_in_text_bv;

    Item<sdsl::int_vector<>> tail_idxs_by_heads_in_text;

    using BV2 = sdsl::sd_vector<>;
//    using BV2 = sdsl::bit_vector;
    Item<BV2> sampled_tails_idx_bv;
    Item<BV2::rank_1_type> rank_sampled_tails_idx_bv;

//    using BV3 = sdsl::sd_vector<>;
    using BV3 = sdsl::bit_vector;
    Item<BV3> marked_sampled_idxs_bv;
    Item<BV3::rank_0_type> rank_marked_sampled_idxs_bv;

    Item<sdsl::int_vector<>> head_marked_sample_trusted_areas;
  };

  std::unordered_map<std::size_t, RIndexComponents> r_index_packs_;

};

#endif //RI_BENCHMARK_FACTORY_H_
