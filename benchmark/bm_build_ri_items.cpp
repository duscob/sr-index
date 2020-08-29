//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/15/20.
//

#include <iostream>
#include <algorithm>

#include <gflags/gflags.h>

#include <benchmark/benchmark.h>

#include <sdsl/config.hpp>
#include <sdsl/construct.hpp>

#include <ri/bwt.h>
#include <ri/rle_string.hpp>
#include "definitions.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_bool(rebuild, false, "Rebuild all the items.");
DEFINE_bool(sais, true, "SE_SAIS or LIBDIVSUFSORT algorithm for Suffix Array construction.");

void SetupCommonCounters(benchmark::State &t_state) {
  t_state.counters["n"] = 0;
  t_state.counters["r"] = 0;
  t_state.counters["s"] = 0;
  t_state.counters["r'"] = 0;
  t_state.counters["mr'"] = 0;
}

//auto BM_Template = [](benchmark::State &t_state, auto *t_config, const auto &_idx) {
//  for (auto _ : t_state) {
//  }
//
//  SetupCommonCounters(t_state);
//};

/// Build text representation to use with SDSL functionalities.
auto BM_BuildText = [](benchmark::State &t_state, auto *t_config, const auto &t_data_path) {
  std::size_t n;
  for (auto _ : t_state) {
    sdsl::int_vector<8> text;
    {
      // Load input text by streaming from disc
      std::string input;
      {
        std::ifstream fs(t_data_path);
        std::stringstream buffer;
        buffer << fs.rdbuf();

        input = buffer.str();
      }

      // Construct text representation for SDSL use.
      n = input.size();
      text.resize(input.size() + 1);

      replace_copy(input.begin(), input.end(), text.begin(), 0, 2);

      text[text.size() - 1] = 0; // Append symbol zero at the end
    }

    sdsl::store_to_cache(text, sdsl::conf::KEY_TEXT, *t_config);
//    sdsl::util::clear(text);
  }

  SetupCommonCounters(t_state);
  t_state.counters["n"] = n;
};

/// Build Suffix Array.
auto BM_BuildSA = [](benchmark::State &t_state, auto *t_config) {
  for (auto _ : t_state) {
    // Use SDSL functionality to build the SA
    sdsl::construct_sa<8>(*t_config);
  }

  SetupCommonCounters(t_state);
};

/// Build Burrows Wheeler Transforms and compute its runs.
auto BM_BuildBWT = [](benchmark::State &t_state, auto *t_config) {
  typedef sdsl::int_vector<>::size_type size_type;
  typedef sdsl::int_vector<8> text_type;
  typedef sdsl::int_vector_buffer<8> bwt_type;

  // Load text from disk
  sdsl::int_vector<8> text;
  load_from_cache(text, sdsl::conf::KEY_TEXT, *t_config);
  auto n = text.size();
  uint8_t bwt_width = text.width();

  sdsl::int_vector<> bwt_heads_pos; // BWT run heads positions in BWT Array
  sdsl::int_vector<> bwt_heads_text_pos; // BWT run heads positions in Text
  std::vector<std::size_t> bwt_heads_text_pos_vec; // BWT run heads positions in Text (vector)
  sdsl::int_vector<> bwt_tails_pos; // BWT run tails positions in BWT Array
  sdsl::int_vector<> bwt_tails_text_pos; // BWT run tails positions in Text
  std::vector<std::size_t> bwt_tails_text_pos_vec; // BWT run tails positions in Text (vector)
  sdsl::bit_vector bwt_heads_in_text_bv; // Mark positions of BWT run heads in text (bitvector)
  sdsl::sd_vector<> bwt_heads_in_text_bv_sd; // Mark positions of BWT run heads in text (sparse bitvector)
  std::vector<std::size_t> f(256, 0); // F Array

  std::size_t r;

  for (auto _ : t_state) {
    // Prepare to stream SA from disc
    size_type buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
    sdsl::int_vector_buffer<> sa_buf(sdsl::cache_file_name(sdsl::conf::KEY_SA, *t_config), std::ios::in, buffer_size);

    // Build BWT sequentially by streaming SA and random access to text
    std::string bwt_file = cache_file_name(sdsl::conf::KEY_BWT, *t_config);
    bwt_type bwt_buf(bwt_file, std::ios::out, buffer_size, bwt_width);

    f.clear();
    f.resize(256, 0);

    auto report_bwt = [&bwt_buf, &f](auto idx, auto symbol) {
      bwt_buf.push_back(symbol);
      ++f[symbol + 1];
    };

    // Report BWT run heads
    std::vector<std::size_t> bwt_heads_pos_vec;
    bwt_heads_text_pos_vec.clear();
    auto report_bwt_head =
        [&bwt_heads_pos_vec, &bwt_heads_text_pos_vec](
            const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
          bwt_heads_pos_vec.emplace_back(idx);
          bwt_heads_text_pos_vec.emplace_back(pos);
        };

    // Report BWT run heads
    std::vector<std::size_t> bwt_tails_pos_vec;
    bwt_tails_text_pos_vec.clear();
    auto report_bwt_tail =
        [&bwt_tails_pos_vec, &bwt_tails_text_pos_vec](
            const auto &run, const auto &idx, const auto &symbol, const auto &pos) {
          bwt_tails_pos_vec.emplace_back(idx);
          bwt_tails_text_pos_vec.emplace_back(pos);
        };

    // Compute BWT and its runs
    r = ri::computeBWT(text.size(),
                       [&sa_buf](auto idx) { return sa_buf[idx]; },
                       [&text](auto idx) { return text[idx]; },
                       report_bwt,
                       report_bwt_head,
                       report_bwt_tail);

    assert(text.size() == bwt_buf.size());
    assert(r == bwt_heads_pos_vec.size());
    assert(r == bwt_tails_pos_vec.size());

    auto log_n = sdsl::bits::hi(n) + 1;

    // Build compact representations for heads and tails of the BWT runs
    bwt_heads_pos = sdsl::int_vector<>(r, 0, log_n);
    bwt_heads_text_pos = sdsl::int_vector<>(r, 0, log_n);

    bwt_tails_pos = sdsl::int_vector<>(r, 0, log_n);
    bwt_tails_text_pos = sdsl::int_vector<>(r, 0, log_n);

    bwt_heads_in_text_bv = sdsl::bit_vector(text.size(), 0);

    for (std::size_t i = 0; i < r; ++i) {
      bwt_heads_pos[i] = bwt_heads_pos_vec[i];
      bwt_heads_text_pos[i] = bwt_heads_text_pos_vec[i];
      bwt_heads_in_text_bv[bwt_heads_text_pos_vec[i]] = 1;

      bwt_tails_pos[i] = bwt_tails_pos_vec[i];
      bwt_tails_text_pos[i] = bwt_tails_text_pos_vec[i];
    }
    bwt_heads_in_text_bv_sd = sdsl::sd_vector<>(bwt_heads_in_text_bv);

    // Build F array.
    for (std::size_t i = 1; i < 256; ++i) {
      f[i] += f[i - 1];
    }

    assert(f[255] == n);
  }

  sdsl::store_to_cache(bwt_heads_pos, ri::KEY_BWT_HEADS, *t_config);
  sdsl::store_to_cache(bwt_heads_text_pos, ri::KEY_BWT_HEADS_TEXT_POS, *t_config);
  sdsl::store_to_cache(bwt_heads_text_pos_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  sdsl::store_to_cache(bwt_heads_in_text_bv, ri::KEY_BWT_HEADS_TEXT_POS + "_bv", *t_config);
  sdsl::store_to_cache(bwt_heads_in_text_bv_sd, ri::KEY_BWT_HEADS_TEXT_POS + "_bv_sd", *t_config);

  sdsl::store_to_cache(bwt_tails_pos, ri::KEY_BWT_TAILS, *t_config);
  sdsl::store_to_cache(bwt_tails_text_pos, ri::KEY_BWT_TAILS_TEXT_POS, *t_config);
  sdsl::store_to_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  sdsl::store_to_cache(f, ri::KEY_F, *t_config);

  SetupCommonCounters(t_state);
  t_state.counters["r"] = r;
};

/// Build run length codification of Burrows Wheeler Transforms.
auto BM_BuildBWTRunLengthEncoded = [](benchmark::State &t_state, auto *t_config) {
  std::size_t buffer_size = 1000000; // buffer_size is a multiple of 8!, TODO: still true?
  std::string bwt_file = cache_file_name(sdsl::conf::KEY_BWT, *t_config);
  sdsl::int_vector_buffer<8> bwt_buf(bwt_file, std::ios::in, buffer_size);

  std::string bwt_s;
  replace_copy(bwt_buf.begin(), bwt_buf.end(), back_inserter(bwt_s), 0, 1);

  ri::rle_string<> bwt_rle;
  for (auto _ : t_state) {
    bwt_rle = ri::rle_string<>(bwt_s);
  }

  sdsl::store_to_cache(bwt_rle, ri::KEY_BWT_RLE, *t_config);

  SetupCommonCounters(t_state);
};

/// Sort BWT run tails by its positions in the text.
auto BM_SortBWTTailsTextPos = [](benchmark::State &t_state, auto *t_config) {
  std::vector<std::size_t> bwt_tails_text_pos_vec;
  sdsl::load_from_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  std::vector<std::size_t> tails_idxs; // Indices of the run tails sorted by its text positions

  for (auto _ : t_state) {
    tails_idxs.clear();
    tails_idxs.resize(bwt_tails_text_pos_vec.size());
    iota(tails_idxs.begin(), tails_idxs.end(), 0);

    sort(tails_idxs.begin(),
         tails_idxs.end(),
         [&bwt_tails_text_pos_vec](const auto &a, const auto &b) -> bool {
           return bwt_tails_text_pos_vec[a] < bwt_tails_text_pos_vec[b];
         });
  }

  sdsl::store_to_cache(tails_idxs, ri::KEY_BWT_TAILS_TEXT_POS_SORTED_IDX + "_vec", *t_config);

  SetupCommonCounters(t_state);
};

/// Sort BWT run heads by its positions in the text.
auto BM_SortBWTHeadsTextPos = [](benchmark::State &t_state, auto *t_config) {
  std::vector<std::size_t> bwt_heads_text_pos_vec;
  sdsl::load_from_cache(bwt_heads_text_pos_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  auto r = bwt_heads_text_pos_vec.size();
  auto log_r = sdsl::bits::hi(r) + 1;

  std::vector<std::size_t> heads_idxs; // Indices of the run tails sorted by its text positions
  sdsl::int_vector<> tail_idxs_by_heads_in_text;

  for (auto _ : t_state) {
    heads_idxs.clear();
    heads_idxs.resize(r);
    iota(heads_idxs.begin(), heads_idxs.end(), 0);

    sort(heads_idxs.begin(),
         heads_idxs.end(),
         [&bwt_heads_text_pos_vec](const auto &a, const auto &b) -> bool {
           return bwt_heads_text_pos_vec[a] < bwt_heads_text_pos_vec[b];
         });

    tail_idxs_by_heads_in_text = sdsl::int_vector<>(r, 0, log_r);
    transform(heads_idxs.begin(),
              heads_idxs.end(),
              tail_idxs_by_heads_in_text.begin(),
              [r](auto i) { return (i + r - 1) % r; });
  }

  sdsl::store_to_cache(heads_idxs, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  sdsl::store_to_cache(tail_idxs_by_heads_in_text,
                       ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT,
                       *t_config);

  SetupCommonCounters(t_state);
};

/// Build BWT run tails sampling.
auto BM_BuildBWTTailsSampling = [](benchmark::State &t_state, auto *t_config) {
  std::size_t s = t_state.range(0); // Sampling size

  std::vector<std::size_t> bwt_tails_text_pos_vec;
  sdsl::load_from_cache(bwt_tails_text_pos_vec, ri::KEY_BWT_TAILS_TEXT_POS + "_vec", *t_config);

  std::vector<std::size_t> tails_idxs_sorted;
  sdsl::load_from_cache(tails_idxs_sorted, ri::KEY_BWT_TAILS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  std::size_t r = tails_idxs_sorted.size();

  // We must sample the run-tails previous to the first and last run-head in the text
  std::size_t prev_tail_to_first_head_in_text;
  std::size_t prev_tail_to_last_head_in_text;
  {
    std::vector<std::size_t> bwt_head_idxs_sorted_in_text;
    sdsl::load_from_cache(bwt_head_idxs_sorted_in_text, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);

    assert(bwt_head_idxs_sorted_in_text.size() == r);

    prev_tail_to_first_head_in_text = (bwt_head_idxs_sorted_in_text.front() + r - 1) % r;
    prev_tail_to_last_head_in_text = (bwt_head_idxs_sorted_in_text.back() + r - 1) % r;
  }

  std::vector<std::size_t> sampled_idxs_vec; // Indices of sampled BWT run tails
  sdsl::int_vector<> sampled_tails_in_text; // Text position of sampled BWT run tails in order of the BWT Array
  sdsl::bit_vector sampled_tails_idx_bv; // Mark which tails are sampled (bitvector)
  sdsl::sd_vector<> sampled_tails_idx_bv_sd; // Mark which tails are sampled (sparse bitvector)

  std::size_t r_prime = 0;

  for (auto _ : t_state) {
    sampled_idxs_vec.clear();
    sampled_idxs_vec.reserve(tails_idxs_sorted.size() / 2);

    // Compute sampled BWT run tails
    auto last_sampled_idx = tails_idxs_sorted[0];
    auto last_sampled = bwt_tails_text_pos_vec[last_sampled_idx];
    sampled_idxs_vec.emplace_back(last_sampled_idx);

    std::size_t prev_idx = tails_idxs_sorted[1];
    std::size_t prev = bwt_tails_text_pos_vec[prev_idx];
    for (std::size_t i = 2; i < tails_idxs_sorted.size(); ++i) {
      auto current_idx = tails_idxs_sorted[i];
      auto current = bwt_tails_text_pos_vec[current_idx];

      if (prev_idx == prev_tail_to_first_head_in_text ||
          prev_idx == prev_tail_to_last_head_in_text ||
          s < current - last_sampled) {
        // Sample the previous run tails in text
        sampled_idxs_vec.emplace_back(prev_idx);
        last_sampled = prev;
      }

      prev_idx = current_idx;
      prev = current;
    }

    sampled_idxs_vec.emplace_back(prev_idx);
    r_prime = sampled_idxs_vec.size();

    // Sort indexes of sampled BWT run tails in order of the BWT Array
    sort(sampled_idxs_vec.begin(), sampled_idxs_vec.end());

    auto log_n = sdsl::bits::hi(prev) + 1;
    sampled_tails_in_text = sdsl::int_vector<>(r_prime, 0, log_n);

    sampled_tails_idx_bv = sdsl::bit_vector(r, 0);

    // Compute text position of the sampled BWT run tails and mark the sampled tails indices in a bitvector
    for (std::size_t i = 0; i < r_prime; ++i) {
      sampled_tails_in_text[i] = bwt_tails_text_pos_vec[sampled_idxs_vec[i]];
      sampled_tails_idx_bv[sampled_idxs_vec[i]] = 1;
    }

    sampled_tails_idx_bv_sd = sdsl::sd_vector<>(sampled_tails_idx_bv);
  }

  sdsl::store_to_cache(sampled_tails_in_text, std::to_string(s) + "_" + ri::KEY_BWT_TAILS_TEXT_POS_SAMPLED, *t_config);
  sdsl::store_to_cache(sampled_idxs_vec, std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_vec", *t_config);
  sdsl::store_to_cache(sampled_tails_idx_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv",
                       *t_config);
  sdsl::store_to_cache(sampled_tails_idx_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_bv_sd",
                       *t_config);

  SetupCommonCounters(t_state);
  t_state.counters["s"] = s;
  t_state.counters["r'"] = r_prime;
};

/// Build BWT run heads sampling.
auto BM_BuildBWTHeadsSampling = [](benchmark::State &t_state, auto *t_config) {
  std::size_t s = t_state.range(0);

  std::vector<std::size_t> heads_in_text_vec;
  sdsl::load_from_cache(heads_in_text_vec, ri::KEY_BWT_HEADS_TEXT_POS + "_vec", *t_config);
  std::size_t r = heads_in_text_vec.size();

  std::vector<std::size_t> sampled_tail_idxs_vec;
  sdsl::load_from_cache(sampled_tail_idxs_vec,
                        std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX + "_vec",
                        *t_config);
  std::size_t r_prime = sampled_tail_idxs_vec.size();
  auto log_r_prime = sdsl::bits::hi(r_prime) + 1;

  std::vector<std::size_t> head_idxs_sorted_in_text;
  sdsl::load_from_cache(head_idxs_sorted_in_text, ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", *t_config);
  assert(head_idxs_sorted_in_text.size() == r);

  // Indices of sampled BWT run tails indices, sorted by the text position of its following run heads (vector)
  std::vector<std::size_t> sampled_tail_idx_by_heads_in_text_vec;
  // Indices of sampled BWT run tails indices, sorted by the text position of its following run heads
  sdsl::int_vector<> sampled_tail_idx_by_heads_in_text;
  // Marked trustful indices of the sampled BWT run heads: bv[i] = true iff between sampled run heads i and i + 1 none run head was deleted
  std::vector<bool> marked_sampled_idxs_vec;
  sdsl::bit_vector marked_sampled_idxs_bv;
  sdsl::sd_vector<> marked_sampled_idxs_bv_sd;
  sdsl::bit_vector sampled_heads_in_text_bv; // Mark positions of sampled BWT run heads in text (bitvector)
  sdsl::sd_vector<> sampled_heads_in_text_bv_sd; // Mark positions of sampled BWT run heads in text (sparse bitvector)

  for (auto _ : t_state) {
    sampled_tail_idx_by_heads_in_text_vec.clear();
    sampled_tail_idx_by_heads_in_text_vec.resize(r_prime);
    std::iota(sampled_tail_idx_by_heads_in_text_vec.begin(), sampled_tail_idx_by_heads_in_text_vec.end(), 0);

    // Sort indices of sampled BWT run tails indices by text position of the next run heads in BWT Array
    sort(sampled_tail_idx_by_heads_in_text_vec.begin(),
         sampled_tail_idx_by_heads_in_text_vec.end(),
         [&heads_in_text_vec, &sampled_tail_idxs_vec, r](const auto &a, const auto &b) {
           auto idx_a = (sampled_tail_idxs_vec[a] + 1) % r;
           auto idx_b = (sampled_tail_idxs_vec[b] + 1) % r;
           return heads_in_text_vec[idx_a] < heads_in_text_vec[idx_b];
         });

    assert((sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec.front()] + 1) % r
               == head_idxs_sorted_in_text.front());
    assert((sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec.back()] + 1) % r
               == head_idxs_sorted_in_text.back());

    marked_sampled_idxs_vec.clear();
    marked_sampled_idxs_vec.resize(r_prime, true);
    marked_sampled_idxs_bv = sdsl::bit_vector(r_prime, true);

    auto n = heads_in_text_vec[head_idxs_sorted_in_text.back()] + 1;
    sampled_heads_in_text_bv = sdsl::bit_vector(n, 0);

    // Marks sampled BWT run heads indices in text and if they are trustworthy
    {
      sampled_heads_in_text_bv[heads_in_text_vec[head_idxs_sorted_in_text.front()]] = true;
      std::size_t prev_idx = 0;
      auto next_sampled_head_idx = (sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec[prev_idx + 1]] + 1) % r;

      for (std::size_t i = 1; i < r - 1; ++i) {
        if (next_sampled_head_idx == head_idxs_sorted_in_text[i]) {
          sampled_heads_in_text_bv[heads_in_text_vec[next_sampled_head_idx]] = true;
          ++prev_idx;
          next_sampled_head_idx = (sampled_tail_idxs_vec[sampled_tail_idx_by_heads_in_text_vec[prev_idx + 1]] + 1) % r;
        } else {
          marked_sampled_idxs_vec[prev_idx] = false;
          marked_sampled_idxs_bv[prev_idx] = false;
        }
      }

      sampled_heads_in_text_bv[heads_in_text_vec[head_idxs_sorted_in_text.back()]] = true;
    }
    sampled_heads_in_text_bv_sd = sdsl::sd_vector<>(sampled_heads_in_text_bv);
    marked_sampled_idxs_bv_sd = sdsl::sd_vector<>(marked_sampled_idxs_bv);

    sampled_tail_idx_by_heads_in_text = sdsl::int_vector<>(r_prime, 0, log_r_prime);
    copy(sampled_tail_idx_by_heads_in_text_vec.begin(),
         sampled_tail_idx_by_heads_in_text_vec.end(),
         sampled_tail_idx_by_heads_in_text.begin());
  }

  sdsl::store_to_cache(sampled_tail_idx_by_heads_in_text,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_SAMPLED_IDX_BY_HEAD_IN_TEXT,
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_vec,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_vec",
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv",
                       *t_config);
  sdsl::store_to_cache(marked_sampled_idxs_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_TAILS_MARKED_SAMPLED_IDX_BY_HEAD_IN_TEXT + "_bv_sd",
                       *t_config);
  sdsl::store_to_cache(sampled_heads_in_text_bv,
                       std::to_string(s) + "_" + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv",
                       *t_config);
  sdsl::store_to_cache(sampled_heads_in_text_bv_sd,
                       std::to_string(s) + "_" + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv_sd",
                       *t_config);

  SetupCommonCounters(t_state);
  t_state.counters["s"] = s;
  t_state.counters["r'"] = r_prime;
  t_state.counters["mr'"] = count(marked_sampled_idxs_vec.begin(), marked_sampled_idxs_vec.end(), true);
};

int main(int argc, char **argv) {
  gflags::SetUsageMessage("This program calculates the ri items for the given text.");
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  if (FLAGS_data.empty()) {
    std::cerr << "Command-line error!!!" << std::endl;
    return 1;
  }

  sdsl::construct_config::byte_algo_sa = FLAGS_sais
                                         ? sdsl::SE_SAIS
                                         : sdsl::LIBDIVSUFSORT; // or LIBDIVSUFSORT for less space-efficient but faster construction

  std::string data_path = FLAGS_data;

  sdsl::cache_config config(false, ".", sdsl::util::basename(FLAGS_data));

  if (!cache_file_exists(sdsl::conf::KEY_TEXT, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildText", BM_BuildText, &config, data_path);
  }

  if (!cache_file_exists(sdsl::conf::KEY_SA, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildSA", BM_BuildSA, &config);
  }

  if (!cache_file_exists(sdsl::conf::KEY_BWT, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildBWT", BM_BuildBWT, &config);
  }

  if (!cache_file_exists(ri::KEY_BWT_RLE, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildBWTRLE", BM_BuildBWTRunLengthEncoded, &config);
  }

  if (!cache_file_exists(ri::KEY_BWT_TAILS_TEXT_POS_SORTED_IDX + "_vec", config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("SortBWTTailsTPos", BM_SortBWTTailsTextPos, &config);
  }

  if (!cache_file_exists(ri::KEY_BWT_HEADS_TEXT_POS_SORTED_IDX + "_vec", config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("SortBWTHeadsTPos", BM_SortBWTHeadsTextPos, &config);
  }

  if (!cache_file_exists("16_" + ri::KEY_BWT_TAILS_TEXT_POS_SAMPLED, config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildBWTTailsSampling", BM_BuildBWTTailsSampling, &config)
        ->RangeMultiplier(2)
        ->Range(4, 2u << 8u);
  }

  if (!cache_file_exists("16_" + ri::KEY_BWT_HEADS_SAMPLED_TEXT_POS + "_bv", config) || FLAGS_rebuild) {
    benchmark::RegisterBenchmark("BuildBWTHeadsSampling", BM_BuildBWTHeadsSampling, &config)
        ->RangeMultiplier(2)
        ->Range(4, 2u << 8u);
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();

  return 0;
}
