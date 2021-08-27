//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/11/21.
//

#ifndef SRI_CONSTRUCT_H_
#define SRI_CONSTRUCT_H_

#include <sdsl/config.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/memory_management.hpp>

#include "psi.h"
#include "tools.h"
#include "phi.h"
#include "rle_string.hpp"
#include "bwt.h"

namespace sri {

template<uint8_t t_width>
struct key_trait {
  static const std::string KEY_BWT;
  static const std::string KEY_BWT_RLE;

  static const std::string KEY_BWT_RUN_FIRST;
  static const std::string KEY_BWT_RUN_FIRST_TEXT_POS;
//  static const std::string KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX;

  static const std::string KEY_BWT_RUN_LAST;
  static const std::string KEY_BWT_RUN_LAST_TEXT_POS;
  static const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX;
  static const std::string KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX;

  static const std::string KEY_ALPHABET;
};

template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT = sdsl::key_bwt_trait<t_width>::KEY_BWT;
template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RLE = key_trait<t_width>::KEY_BWT + "_rle";

template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_FIRST = key_trait<t_width>::KEY_BWT + "_run_first";
template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS = key_trait<t_width>::KEY_BWT_RUN_FIRST + "_text_pos";

template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_LAST = key_trait<t_width>::KEY_BWT + "_run_last";
template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS = key_trait<t_width>::KEY_BWT_RUN_LAST + "_text_pos";
template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX =
    key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_idx";
template<uint8_t t_width>
const std::string key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX =
    key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS + "_sorted_to_first_idx";

template<>
const std::string key_trait<8>::KEY_ALPHABET = "alphabet";
template<>
const std::string key_trait<0>::KEY_ALPHABET = "alphabet_int";

template<uint8_t t_width>
struct alphabet_trait {
  typedef sdsl::byte_alphabet type;
};

template<>
struct alphabet_trait<0> {
  typedef sdsl::int_alphabet<> type;
};

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

  store_to_cache(text, KEY_TEXT, t_config);
}

template<uint8_t t_width>
void constructBWTRLE(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructBWTRLE: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));

  std::string bwt_s;
  replace_copy(bwt_buf.begin(), bwt_buf.end(), back_inserter(bwt_s), 0, 1);

  sri::rle_string<> bwt_rle(bwt_s);

  sdsl::store_to_cache(bwt_rle, key_trait<t_width>::KEY_BWT_RLE, t_config);
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
  auto out_int_vector_buffer = [buffer_size, n_width](const auto &tt_key) {
    return sdsl::int_vector_buffer<>(tt_key, std::ios::out, buffer_size, n_width);
  };

  auto file_bwt_run_first = cache_file_name(key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);
  auto file_bwt_run_first_text_pos = cache_file_name(key_trait<t_width>::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);
  auto file_bwt_run_last = cache_file_name(key_trait<t_width>::KEY_BWT_RUN_LAST, t_config);
  auto file_bwt_run_last_text_pos = cache_file_name(key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_config);

  auto bwt_run_first_pos = out_int_vector_buffer(file_bwt_run_first); // BWT run heads positions in BWT array
  auto bwt_run_first_text_pos = out_int_vector_buffer(file_bwt_run_first_text_pos); // BWT run heads positions in text
  auto bwt_run_last_pos = out_int_vector_buffer(file_bwt_run_last); // BWT run tails positions in BWT array
  auto bwt_run_last_text_pos = out_int_vector_buffer(file_bwt_run_last_text_pos); // BWT run tails positions in text

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
  register_cache_file(file_bwt_run_first, t_config);

  bwt_run_first_text_pos.close();
  register_cache_file(file_bwt_run_first_text_pos, t_config);

  bwt_run_last_pos.close();
  register_cache_file(file_bwt_run_last, t_config);

  bwt_run_last_text_pos.close();
  register_cache_file(file_bwt_run_last_text_pos, t_config);
}

template<uint8_t t_width>
void constructAlphabet(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructAlphabet: width must be `0` for integer alphabet and `8` for byte alphabet");

  sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));
  auto n = bwt_buf.size();

  typename alphabet_trait<t_width>::type alphabet(bwt_buf, n);

  store_to_cache(alphabet, key_trait<t_width>::KEY_ALPHABET, t_config);
}

template<uint8_t t_width>
void constructPsi(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructPsi: width must be `0` for integer alphabet and `8` for byte alphabet");

  typename alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, key_trait<t_width>::KEY_ALPHABET, t_config);

  sdsl::int_vector<> psi;
  {
    sdsl::int_vector_buffer<t_width> bwt_buf(sdsl::cache_file_name(sdsl::key_bwt_trait<t_width>::KEY_BWT, t_config));

    psi = constructPsi(bwt_buf, alphabet);
    sdsl::util::bit_compress(psi);
    sdsl::store_to_cache(psi, sdsl::conf::KEY_PSI, t_config);

    // TODO Don't create psi_enc
    {
      sdsl::enc_vector<> psi_enc(psi);
      sdsl::store_to_cache(psi_enc, sdsl::conf::KEY_PSI, t_config, true);
    }
  }

  // TODO Don't create unused PsiCore
  {
    sri::PsiCoreBV<> psi_core(alphabet.C, psi);
    sdsl::store_to_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  }
  {
    sri::PsiCoreBV<sdsl::sd_vector<>> psi_core(alphabet.C, psi);
    sdsl::store_to_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  }
  {
    sri::PsiCoreBV<sdsl::rrr_vector<>> psi_core(alphabet.C, psi);
    sdsl::store_to_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  }

  {
    sri::PsiCoreRLE<> psi_rle(alphabet.C, psi);
    sdsl::store_to_cache(psi_rle, sdsl::conf::KEY_PSI, t_config, true);
  }
}

template<typename TRAContainer, typename TGetLink>
auto constructMarkToSampleLinks(const TRAContainer &t_mark_text_pos, const TGetLink &t_get_link) {
  auto r = t_mark_text_pos.size();
  auto log_r = sdsl::bits::hi(r) + 1;

  sdsl::int_vector<> mark_idxs(r, 0, log_r); // Indices of the mark values sorted by its text positions
  iota(mark_idxs.begin(), mark_idxs.end(), 0);

  sort(mark_idxs.begin(),
       mark_idxs.end(),
       [&t_mark_text_pos](const auto &a, const auto &b) -> bool {
         return t_mark_text_pos[a] < t_mark_text_pos[b];
       });

  sdsl::int_vector<> mark_to_sample_link(r, 0, log_r);
  transform(mark_idxs.begin(), mark_idxs.end(), mark_to_sample_link.begin(), t_get_link);

  return std::make_pair(std::move(mark_idxs), std::move(mark_to_sample_link));
}

template<uint8_t t_width>
void constructMarkToSampleLinksForPhiForward(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructMarkToSampleLinksForPhiForward: width must be `0` for integer alphabet and `8` for byte alphabet");

  // Marks
  sdsl::int_vector<> bwt_run_last_text_pos; // BWT run tails positions in text
  sdsl::load_from_cache(bwt_run_last_text_pos, key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS, t_config);

  // Mark positions
  sdsl::int_vector<> bwt_run_last;
  sdsl::load_from_cache(bwt_run_last, key_trait<t_width>::KEY_BWT_RUN_LAST, t_config);

  // LF
  sri::rle_string<> bwt_rle;
  sdsl::load_from_cache(bwt_rle, key_trait<t_width>::KEY_BWT_RLE, t_config);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  typename alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, key_trait<t_width>::KEY_ALPHABET, t_config);
  auto n = alphabet.C[alphabet.sigma];

  auto get_f = [&alphabet](auto tt_symbol) { return alphabet.C[alphabet.char2comp[tt_symbol]]; };
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  // Psi
  sri::PsiCoreBV<sdsl::sd_vector<>> psi_core;
  sdsl::load_from_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };

  auto get_c = [&alphabet](auto tt_index) { return sri::computeCForSAIndex(alphabet.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  // Samples
  sdsl::int_vector<> bwt_run_first; // BWT run heads positions in text
  sdsl::load_from_cache(bwt_run_first, key_trait<t_width>::KEY_BWT_RUN_FIRST, t_config);

  auto rank_sample = [&bwt_run_first](const auto &tt_k) {
    return std::lower_bound(bwt_run_first.begin(), bwt_run_first.end(), tt_k) - bwt_run_first.begin();
  };

  auto get_link = [&bwt_run_last, &n, &lf, &psi, &rank_sample](const auto &tt_mark_idx) {
    return sri::computeMarkToSampleLinkForPhiForward(bwt_run_last[tt_mark_idx], n, lf, psi, rank_sample);
  };

  // Compute links
  auto[sorted_mark_idxs, mark_to_sample_links] = constructMarkToSampleLinks(bwt_run_last_text_pos, get_link);

  sdsl::store_to_cache(sorted_mark_idxs, key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX, t_config);
  sdsl::store_to_cache(mark_to_sample_links,
                       key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX,
                       t_config);
}

template<uint8_t t_width, typename TIndex>
void constructSRI(TIndex &t_index, const std::string &t_data_path, sdsl::cache_config &t_config) {

  {
    // Parse Text
    auto event = sdsl::memory_monitor::event("Text");
    const char *KEY_TEXT = sdsl::key_text_trait<t_width>::KEY_TEXT;
    if (!cache_file_exists(KEY_TEXT, t_config)) {
      sri::constructText<8>(t_data_path, t_config);
    }
  }

  {
    // Construct Suffix Array
    auto event = sdsl::memory_monitor::event("SA");
    if (!cache_file_exists(sdsl::conf::KEY_SA, t_config)) {
      sdsl::construct_sa<t_width>(t_config);
    }
  }

  {
    // Construct BWT
    auto event = sdsl::memory_monitor::event("BWT");
    if (!cache_file_exists(key_trait<t_width>::KEY_BWT, t_config)) {
      sdsl::construct_bwt<t_width>(t_config);
    }
  }

  {
    // Construct BWT RLE
    auto event = sdsl::memory_monitor::event("BWT RLE");
    if (!cache_file_exists(key_trait<t_width>::KEY_BWT_RLE, t_config)) {
      constructBWTRLE<t_width>(t_config);
    }
  }

  {
    // Construct BWT Runs
    auto event = sdsl::memory_monitor::event("BWT Runs");
    const auto KEY_BWT_RUN_FIRST = key_trait<t_width>::KEY_BWT_RUN_FIRST;
    if (!cache_file_exists(KEY_BWT_RUN_FIRST, t_config)) {
      constructBWTRuns<t_width>(t_config);
    }
  }

  {
    // Construct Alphabet
    auto event = sdsl::memory_monitor::event("Alphabet");
    if (!cache_file_exists(key_trait<t_width>::KEY_ALPHABET, t_config)) {
      constructAlphabet<t_width>(t_config);
    }
  }

  {
    // Construct Psi
    auto event = sdsl::memory_monitor::event("Psi");
    if (!cache_file_exists(sdsl::conf::KEY_PSI, t_config)) {
      constructPsi<t_width>(t_config);
    }
  }

  {
    // Construct Links from Mark to Sample
    auto event = sdsl::memory_monitor::event("Mark2Sample Links");
    if (!cache_file_exists(key_trait<t_width>::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config)) {
      constructMarkToSampleLinksForPhiForward<t_width>(t_config);
    }
  }

  construct(t_index, t_config);
}

template<typename TBitVector>
void constructBitVectorFromIntVector(const std::string &t_key, sdsl::cache_config &t_config, std::size_t t_bv_size) {
  sdsl::bit_vector bv_tmp(t_bv_size, 0);

  sdsl::int_vector_buffer<> int_buf(sdsl::cache_file_name(t_key, t_config));
  for (int i = 0; i < int_buf.size(); ++i) {
    bv_tmp[int_buf[i]] = true;
  }

  TBitVector bv(std::move(bv_tmp));
  sdsl::store_to_cache(bv, t_key, t_config, true);

  typename TBitVector::rank_1_type bv_rank(&bv);
  sdsl::store_to_cache(bv_rank, t_key, t_config, true);

  typename TBitVector::select_1_type bv_select(&bv);
  sdsl::store_to_cache(bv_select, t_key, t_config, true);
}

}

#endif //SRI_CONSTRUCT_H_
