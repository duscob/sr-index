//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/7/2023.
//

#ifndef SRI_CONSTRUCT_BASE_H_
#define SRI_CONSTRUCT_BASE_H_

#include <string>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/io.hpp>
#include <sdsl/iterators.hpp>

#include "alphabet.h"
#include "rle_string.hpp"

namespace sri {

namespace conf {
const std::string KEY_ALPHABET = "alphabet";

const std::string KEY_BWT_RLE = "bwt_rle";

//! BWT run heads
const std::string KEY_BWT_RUN_FIRST = "bwt_run_first"; // Positions in the BWT array
const std::string KEY_BWT_RUN_FIRST_IDX = KEY_BWT_RUN_FIRST + "_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS = KEY_BWT_RUN_FIRST + "_text_pos";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_BY_LAST = KEY_BWT_RUN_FIRST_TEXT_POS + "_by_last";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_TO_LAST_IDX = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_to_last_idx";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_VALID_MARK = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_valid_mark";
const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_VALID_AREA = KEY_BWT_RUN_FIRST_TEXT_POS + "_sorted_valid_area";

//! BWT run tails
const std::string KEY_BWT_RUN_LAST = "bwt_run_last"; // Positions in the BWT array
const std::string KEY_BWT_RUN_LAST_IDX = KEY_BWT_RUN_LAST + "_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS = KEY_BWT_RUN_LAST + "_text_pos";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_BY_FIRST = KEY_BWT_RUN_LAST_TEXT_POS + "_by_first";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_to_first_idx";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_MARK = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_valid_mark";
const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_VALID_AREA = KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_valid_area";

const std::string KEY_BWT_RUN_CUMULATIVE_COUNT = "bwt_run_cumulative_count";

//! Psi run heads
const std::string KEY_PSI_RUN_FIRST = "psi_run_first"; // Positions in the Psi array
const std::string KEY_PSI_RUN_FIRST_TEXT_POS = KEY_PSI_RUN_FIRST + "_text_pos";
}

template<uint8_t t_width>
void constructAlphabet(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructAlphabet: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
  auto n = bwt_buf.size();

  typename alphabet_trait<t_width>::type alphabet(bwt_buf, n);

  sdsl::store_to_cache(alphabet, conf::KEY_ALPHABET, t_config);
}

template<uint8_t t_width>
void constructBWTRLE(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructBWTRLE: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));

  {
    typename alphabet_trait<t_width>::type alphabet;
    sdsl::load_from_cache(alphabet, conf::KEY_ALPHABET, t_config);

    auto get_symbol = [&bwt_buf, &alphabet](auto tt_i) { return alphabet.char2comp[bwt_buf[tt_i]]; };

    auto bwt_s = sdsl::random_access_container(get_symbol, bwt_buf.size());

    RLEString<> bwt_rle(bwt_s.begin(), bwt_s.end());

    sdsl::store_to_cache(bwt_rle, conf::KEY_BWT_RLE, t_config);
  }
}

}

#endif //SRI_CONSTRUCT_BASE_H_
