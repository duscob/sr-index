//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/7/2023.
//

#ifndef SRI_CONSTRUCT_BASE_H_
#define SRI_CONSTRUCT_BASE_H_

#include <string>

namespace sri {

namespace conf {
const std::string KEY_ALPHABET = "alphabet";

const std::string KEY_BWT_RLE = "bwt_rle";

const std::string KEY_BWT_RUN_FIRST = "bwt_run_first";
const std::string KEY_BWT_RUN_FIRST_IDX = KEY_BWT_RUN_FIRST + "_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS = KEY_BWT_RUN_FIRST + "_text_pos";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_BY_LAST = KEY_BWT_RUN_FIRST_TEXT_POS + "_by_last";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_TO_LAST_IDX = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_to_last_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_VALID_MARK = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_valid_mark";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_VALID_AREA = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_valid_area";

const std::string KEY_BWT_RUN_LAST = "bwt_run_last";
const std::string KEY_BWT_RUN_LAST_IDX = KEY_BWT_RUN_LAST + "_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS = KEY_BWT_RUN_LAST + "_text_pos";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST = KEY_BWT_RUN_LAST_TEXT_POS + "_by_first";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_to_first_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_valid_mark";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_valid_area";

const std::string KEY_BWT_RUN_CUMULATIVE_COUNT = "bwt_run_cumulative_count";
}

}

#endif //SRI_CONSTRUCT_BASE_H_
