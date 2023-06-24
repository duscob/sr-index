//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/24/2023.
//

#ifndef SRI_CONSTRUCT_SDSL_H_
#define SRI_CONSTRUCT_SDSL_H_

#include <cstdint>
#include <string>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/vectors.hpp>

#include "construct_base.h"

namespace sri {
namespace inner_sdsl {

template<uint8_t t_width>
void constructText(const std::string &t_file, sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructText: width must be `0` for integer alphabet and `8` for byte alphabet");

  using TText = sdsl::int_vector<t_width>;

  const auto KEY_TEXT = sdsl::key_text_trait<t_width>::KEY_TEXT;

  TText text;
  auto num_bytes = t_width / 8;
  load_vector_from_file(text, t_file, num_bytes);

  auto it_zero = std::find(text.begin(), text.end(), (uint64_t) 0);
  if (it_zero == text.end()) {
    sdsl::append_zero_symbol(text);
  } else if (it_zero != text.end() - 1) {
    throw std::logic_error(std::string("Error: File \"") + t_file + "\" contains inner zero symbol.");
  }

  sdsl::store_to_cache(text, KEY_TEXT, t_config);
}

template<uint8_t t_width>
void constructBWTRuns(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructBWTRuns: width must be `0` for integer alphabet and `8` for byte alphabet");

  // Prepare to stream BWT and SA from disc
  // TODO Use int_vector_buffer instead int_vector to process big files
//  sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
  sdsl::int_vector<t_width> bwt_buf;
  sdsl::load_from_cache(bwt_buf, sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config);
//  sdsl::int_vector_buffer<> sa_buf(sdsl::cache_file_name(sdsl::conf::KEY_SA, t_config));
  sdsl::int_vector<> sa_buf;
  sdsl::load_from_cache(sa_buf, sdsl::conf::KEY_SA, t_config);

  const auto n = bwt_buf.size();

  auto get_bwt_text_pos = [&sa_buf, n](auto tt_bwt_idx) {
    auto next_pos = sa_buf[tt_bwt_idx];
    return 0 < next_pos ? next_pos - 1 : n - 1;
  };

  // Prepare to BWT runs to disc
  const std::size_t buffer_size = 1 << 20;
  const std::size_t n_width = sdsl::bits::hi(n) + 1;
  auto out_int_vector_buf = [buffer_size, n_width, &t_config](const auto &tt_key) {
    return sdsl::int_vector_buffer<>(cache_file_name(tt_key, t_config), std::ios::out, buffer_size, n_width);
  };

  auto bwt_run_first_pos = out_int_vector_buf(conf::KEY_BWT_RUN_FIRST); // BWT run head positions in BWT array
  auto bwt_run_first_text_pos = out_int_vector_buf(conf::KEY_BWT_RUN_FIRST_TEXT_POS); // BWT run head positions in text
  auto bwt_run_last_pos = out_int_vector_buf(conf::KEY_BWT_RUN_LAST); // BWT run tail positions in BWT array
  auto bwt_run_last_text_pos = out_int_vector_buf(conf::KEY_BWT_RUN_LAST_TEXT_POS); // BWT run tail positions in text

  size_t n_runs = 0; // # BWT runs

  // First BWT value
  auto bwt_symbol = bwt_buf[0];
  auto text_pos = get_bwt_text_pos(0);

  // First position starts the first BWT run.
  bwt_run_first_pos.push_back(0);
  bwt_run_first_text_pos.push_back(text_pos);

  auto prev_bwt_symbol = bwt_symbol;
  auto prev_text_pos = text_pos;
  for (int i = 1; i < n; ++i) {
    bwt_symbol = bwt_buf[i];
    text_pos = get_bwt_text_pos(i);

    if (bwt_symbol != prev_bwt_symbol) {
      // Last position of the previous BWT run
      bwt_run_last_pos.push_back(i - 1);
      bwt_run_last_text_pos.push_back(prev_text_pos);

      ++n_runs;

      // First position of the current BWT run
      bwt_run_first_pos.push_back(i);
      bwt_run_first_text_pos.push_back(text_pos);

      prev_bwt_symbol = bwt_symbol;
    }

    prev_text_pos = text_pos;
  }

  // Last position ends the last BWT run
  bwt_run_last_pos.push_back(n - 1);
  bwt_run_last_text_pos.push_back(text_pos);

  bwt_run_first_pos.close();
  register_cache_file(conf::KEY_BWT_RUN_FIRST, t_config);

  bwt_run_first_text_pos.close();
  register_cache_file(conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);

  bwt_run_last_pos.close();
  register_cache_file(conf::KEY_BWT_RUN_LAST, t_config);

  bwt_run_last_text_pos.close();
  register_cache_file(conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config);
}

template<uint8_t t_width>
void constructIndexBaseItems(const std::string &t_data_path, sdsl::cache_config &t_config) {
  // Parse Text
  const char *KEY_TEXT = sdsl::key_text_trait<t_width>::KEY_TEXT;
  if (!cache_file_exists(KEY_TEXT, t_config)) {
    auto event = sdsl::memory_monitor::event("Text");
    constructText<t_width>(t_data_path, t_config);
  }

  // Construct Suffix Array
  if (!cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
    auto event = sdsl::memory_monitor::event("SA");
    sdsl::construct_sa<t_width>(t_config);
  }

  // Construct BWT
  if (!cache_file_exists(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config)) {
    auto event = sdsl::memory_monitor::event("BWT");
    sdsl::construct_bwt<t_width>(t_config);
  }

  // Construct BWT Runs
  if (!cache_file_exists(conf::KEY_BWT_RUN_FIRST, t_config)) {
    auto event = sdsl::memory_monitor::event("BWT Runs");
    constructBWTRuns<t_width>(t_config);
  }

  // Construct Alphabet
  if (!cache_file_exists(conf::KEY_ALPHABET, t_config)) {
    auto event = sdsl::memory_monitor::event("Alphabet");
    constructAlphabet<t_width>(t_config);
  }

  // Construct BWT RLE
  if (!cache_file_exists(conf::KEY_BWT_RLE, t_config)) {
    auto event = sdsl::memory_monitor::event("BWT RLE");
    constructBWTRLE<t_width>(t_config);
  }
}

} // namespace inner_sdsl
} // namespace sri

#endif //SRI_CONSTRUCT_SDSL_H_
