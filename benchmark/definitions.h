//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/16/20.
//

#ifndef RI_BENCHMARK_DEFINITIONS_H_
#define RI_BENCHMARK_DEFINITIONS_H_

#include <string>

namespace ri {

const std::string KEY_BWT_HEADS = "bwt_heads";
const std::string KEY_BWT_HEADS_TEXT_POS = "bwt_heads_tpos";
const std::string KEY_BWT_TAILS = "bwt_tails";
const std::string KEY_BWT_TAILS_TEXT_POS = "bwt_tails_tpos";

const std::string KEY_BWT_HEADS_SAMPLED_TEXT_POS = "bwt_heads_sampled_tpos";
const std::string KEY_BWT_TAILS_TEXT_POS_SAMPLED = "bwt_tails_sampled_tpos";
//const std::string KEY_BWT_TAILS_MARKED_SAMPLED_IDX = "bwt_tails_marked_sampled_idx";

const std::string KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT = "bwt_tails_sampled_idx_by_heads_tpos";
const std::string KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT = "bwt_tails_marked_sampled_idx_by_heads_tpos";

const std::string KEY_BWT_HEADS_TEXT_POS_SORTED_IDX = "bwt_heads_tpos_sorted_idx";
const std::string KEY_BWT_TAILS_TEXT_POS_SORTED_IDX = "bwt_tails_tpos_sorted_idx";
const std::string KEY_BWT_TAILS_SAMPLED_IDX = "bwt_tails_sampled_idx";

const std::string KEY_F = "f";

}

#endif //RI_BENCHMARK_DEFINITIONS_H_
