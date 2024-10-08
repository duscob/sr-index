//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/7/2023.
//

#ifndef SRI_CONSTRUCT_BIG_BWT_H_
#define SRI_CONSTRUCT_BIG_BWT_H_

#ifndef BIGBWT_EXE
#define BIGBWT_EXE "bigbwt"
#endif

#include <cstdint>
#include <string>
#include <filesystem>

#include <sdsl/config.hpp>
#include <sdsl/int_vector_buffer.hpp>

#include "construct_base.h"

namespace sri {

namespace conf {
const std::string KEY_BIG_BWT = "big_bwt";
const std::string KEY_BIG_BWT_LOG = "big_bwt_log";
const std::string KEY_BIG_BWT_SSA = "big_bwt_ssa";
const std::string KEY_BIG_BWT_ESA = "big_bwt_esa";
} // namespace conf

namespace inner_big_bwt {

void registerBigBWTFiles(const std::string &t_data_path, sdsl::cache_config &t_config) {
  t_config.file_map[conf::KEY_BIG_BWT] = t_data_path + ".bwt";
  t_config.file_map[conf::KEY_BIG_BWT_LOG] = t_data_path + ".log";
  t_config.file_map[conf::KEY_BIG_BWT_SSA] = t_data_path + ".ssa";
  t_config.file_map[conf::KEY_BIG_BWT_ESA] = t_data_path + ".esa";
}

void constructBWT(const std::string &t_data_path, sdsl::cache_config &t_config, const std::string &t_bigbwt_exe) {
  auto command = t_bigbwt_exe + " -s -e " + t_data_path + " > /dev/null 2>&1";
  std::system(command.c_str());

  // Prepare to stream BWT symbols from disc
  std::ifstream bwt_file(cache_file_name(conf::KEY_BIG_BWT, t_config));
  // Prepare to store BWT and runs to disc
  sdsl::int_vector_buffer<8> bwt_buf(cache_file_name(sdsl::conf::KEY_BWT, t_config), std::ios::out);
  uint8_t symbol;
  while (bwt_file >> symbol) {
    bwt_buf.push_back(symbol);
  }
  register_cache_file(sdsl::conf::KEY_BWT, t_config);
}

void constructBWTRuns(sdsl::cache_config &t_config) {
  const auto n = std::filesystem::file_size(sdsl::cache_file_name(conf::KEY_BIG_BWT, t_config));

  // Prepare to store BWT runs to disc
  const std::size_t buffer_size = 1 << 20;
  const std::size_t n_width = sdsl::bits::hi(n) + 1;
  auto out_int_vector_buf = [buffer_size, n_width, &t_config](const auto &tt_key) {
    return sdsl::int_vector_buffer<>(cache_file_name(tt_key, t_config), std::ios::out, buffer_size, n_width);
  };

  auto read_runs = [n, &t_config, &out_int_vector_buf](
      const auto &tt_key, const auto &tt_key_bwt_run_pos, const auto &tt_key_bwt_run_text_pos
  ) {
    // Prepare to stream BWT run positions <j, SA[j]> from disc
    std::ifstream input(cache_file_name(tt_key, t_config));
    if (!input) return;

    auto bwt_run_pos = out_int_vector_buf(tt_key_bwt_run_pos); // BWT run positions in BWT array
    auto bwt_run_text_pos = out_int_vector_buf(tt_key_bwt_run_text_pos); // BWT run positions in text

    uint64_t j = 0;
    uint64_t sa_j = 0;
    while (input.read((char *) &j, 5) && input.read((char *) &sa_j, 5)) {
      bwt_run_pos.push_back(j);
      bwt_run_text_pos.push_back(sa_j ? sa_j - 1 : n - 1);
    }

    bwt_run_pos.close();
    register_cache_file(tt_key_bwt_run_pos, t_config);

    bwt_run_text_pos.close();
    register_cache_file(tt_key_bwt_run_text_pos, t_config);
  };

  read_runs(conf::KEY_BIG_BWT_SSA, conf::KEY_BWT_RUN_FIRST, conf::KEY_BWT_RUN_FIRST_TEXT_POS);
  read_runs(conf::KEY_BIG_BWT_ESA, conf::KEY_BWT_RUN_LAST, conf::KEY_BWT_RUN_LAST_TEXT_POS);
}

template<uint8_t t_width>
void constructIndexBaseItems(const std::string &t_data_path,
                             sdsl::cache_config &t_config,
                             const std::string &t_big_bwt_exe = BIGBWT_EXE) {

  // Register files generated by big-bwt
  registerBigBWTFiles(t_data_path, t_config);

  // Construct BWT
  if (!cache_file_exists(conf::KEY_BIG_BWT, t_config)) {
    auto event = sdsl::memory_monitor::event("BWT");
    constructBWT(t_data_path, t_config, t_big_bwt_exe);
  }

  // Construct BWT Runs
  if (!cache_file_exists(conf::KEY_BWT_RUN_FIRST, t_config)) {
    auto event = sdsl::memory_monitor::event("BWT Runs");
    constructBWTRuns(t_config);
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

} // namespace inner_big_bwt
} // namespace sri

#endif //SRI_CONSTRUCT_BIG_BWT_H_
