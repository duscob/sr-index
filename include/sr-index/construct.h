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

#include "config.h"
#include "construct_base.h"
#include "construct_sdsl.h"
#include "construct_big_bwt.h"
#include "alphabet.h"
#include "psi.h"
#include "tools.h"
#include "phi.h"
#include "rle_string.hpp"
#include "bwt.h"
#include "io.h"

namespace sri {

template<typename T>
auto construct(const sdsl::int_vector<> &t_iv) {
  return T(t_iv);
}

template<typename TIndex>
void construct(TIndex& t_index, const std::string& t_data_path, Config& t_config) {
  constructItems(t_index, t_config);

  t_index.load(t_config);
}

template<uint8_t t_width>
void constructIndexBaseItems(const std::string &t_data_path, sri::Config &t_config) {
  switch (t_config.sa_algo) {
    case SDSL_LIBDIVSUFSORT:
      sdsl::construct_config().byte_algo_sa = sdsl::LIBDIVSUFSORT;
      inner_sdsl::constructIndexBaseItems<t_width>(t_data_path, t_config);
      break;
    case SDSL_SE_SAIS:
      sdsl::construct_config().byte_algo_sa = sdsl::SE_SAIS;
      inner_sdsl::constructIndexBaseItems<t_width>(t_data_path, t_config);
      break;
    case BIG_BWT:
      inner_big_bwt::constructIndexBaseItems<t_width>(t_data_path, t_config);
      break;
  }
}

inline auto KeySortedByAlphabet(const std::string& t_key) {
  return t_key + "_sorted_alphabet";
}

inline auto getExtremes(const sdsl::cache_config& t_config, const std::string& t_key, bool t_add_type_hash = false) {
  const auto filename = t_add_type_hash
                          ? sdsl::cache_file_name<sdsl::int_vector<>>(t_key, t_config)
                          : sdsl::cache_file_name(t_key, t_config);
  sdsl::int_vector_buffer<> buffer(filename);
  return std::array<std::size_t, 2>{buffer[0], buffer[buffer.size() - 1]};
}

template<uint8_t t_width>
void constructPsi(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructPsi: width must be `0` for integer alphabet and `8` for byte alphabet");

  typename alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, conf::KEY_ALPHABET, t_config);

  sdsl::int_vector<> psi;
  {
    RLEString<> bwt_rle;
    sdsl::load_from_cache(bwt_rle, conf::KEY_BWT_RLE, t_config);
    auto get_bwt_symbol = [&bwt_rle](size_t tt_i) { return bwt_rle[tt_i]; };

    psi = constructPsi(get_bwt_symbol, alphabet.C);
    sdsl::util::bit_compress(psi);
    sdsl::store_to_cache(psi, sdsl::conf::KEY_PSI, t_config);
  }

  {
    sri::PsiCoreRLE<> psi_rle(alphabet.C, psi);
    sri::store_to_cache(psi_rle, sdsl::conf::KEY_PSI, t_config, true);
  }
}

template<uint8_t t_width>
void constructPsiRuns(Config &t_config) {
  static_assert(t_width == 0 || t_width == 8,
                "constructPsiRuns: width must be `0` for integer alphabet and `8` for byte alphabet");

  using namespace conf;

  typename alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, conf::KEY_ALPHABET, t_config);

  RLEString<> bwt_rle;
  sdsl::load_from_cache(bwt_rle, conf::KEY_BWT_RLE, t_config);

  for (const auto &part : {conf::kHead, conf::kTail}) {
    const auto &key_bwt_run_pos = t_config.keys[kBWT][part][kPos];
    const auto &key_bwt_run_text_pos = t_config.keys[kBWT][part][kTextPos];
    const auto &key_psi_run_text_pos = t_config.keys[kPsi][part][kTextPos];

    std::vector<std::vector<std::size_t>> psi_run_text_pos_partial(alphabet.sigma);

    auto bwt_run_pos = sdsl::int_vector_buffer<>(cache_file_name(key_bwt_run_pos, t_config));
    auto bwt_run_text_pos = sdsl::int_vector_buffer<>(cache_file_name(key_bwt_run_text_pos, t_config));

    auto r = bwt_run_pos.size();
    for (std::size_t i = 0; i < r; ++i) {
      psi_run_text_pos_partial[bwt_rle[bwt_run_pos[i]]].emplace_back(bwt_run_text_pos[i]);
    }

    auto psi_run_text_pos = sdsl::int_vector_buffer<>(
      sdsl::cache_file_name<sdsl::int_vector<>>(key_psi_run_text_pos, t_config),
      std::ios::out,
      1 << 20,
      bwt_run_text_pos.width()
    );

    for (const auto &symbol_text_pos : psi_run_text_pos_partial) {
      for (const auto &text_pos : symbol_text_pos) {
        psi_run_text_pos.push_back(text_pos);
      }
    }

    psi_run_text_pos.close();
    sdsl::register_cache_file(key_psi_run_text_pos, t_config);
  }
}

template<typename TRAContainer>
auto sortIndices(const TRAContainer &t_values) {
  auto n = t_values.size();
  auto log_n = sdsl::bits::hi(n) + 1;

  sdsl::int_vector<> values_idx(n, 0, log_n); // Indices of the values sorted
  std::iota(values_idx.begin(), values_idx.end(), 0);

  std::sort(values_idx.begin(),
            values_idx.end(),
            [&t_values](const auto &a, const auto &b) -> bool { return t_values[a] < t_values[b]; });

  return values_idx;
}

template<typename TRAContainer, typename TGetLink>
auto constructMarkToSampleLinks(const TRAContainer &t_marks_text_pos, const TGetLink &t_get_link) {
  // Indices of the mark values sorted by its text positions
  sdsl::int_vector<> marks_idx = sortIndices(t_marks_text_pos);

  auto r = t_marks_text_pos.size();
  auto log_r = sdsl::bits::hi(r) + 1;
  sdsl::int_vector<> mark_to_sample_link(r, 0, log_r);
  std::transform(marks_idx.begin(), marks_idx.end(), mark_to_sample_link.begin(), t_get_link);

  return std::make_pair(std::move(marks_idx), std::move(mark_to_sample_link));
}

template<uint8_t t_width>
void constructMarkToSampleLinksForPhiForwardWithBWTRuns(sdsl::cache_config &t_config) {
  static_assert(t_width == 0 or t_width == 8,
                "constructMarkToSampleLinksForPhiForwardWithBWTRuns: width must be `0` for integer alphabet and `8` "
                "for byte alphabet");

  // Marks
  sdsl::int_vector<> bwt_run_last_text_pos; // BWT run tails positions in text
  sdsl::load_from_cache(bwt_run_last_text_pos, conf::KEY_BWT_RUN_LAST_TEXT_POS, t_config);

  // Mark positions
  sdsl::int_vector<> bwt_run_last;
  sdsl::load_from_cache(bwt_run_last, conf::KEY_BWT_RUN_LAST, t_config);

  // LF
  RLEString<> bwt_rle;
  sdsl::load_from_cache(bwt_rle, conf::KEY_BWT_RLE, t_config);
  auto get_char = sri::buildRandomAccessForContainer(std::cref(bwt_rle));
  auto get_rank_of_char = sri::buildRankOfChar(std::cref(bwt_rle));

  typename alphabet_trait<t_width>::type alphabet;
  sdsl::load_from_cache(alphabet, conf::KEY_ALPHABET, t_config);
  auto n = alphabet.C[alphabet.sigma];

  auto get_f = [&alphabet](auto tt_symbol) { return alphabet.C[tt_symbol]; };
  auto lf = sri::buildBasicLF(get_char, get_rank_of_char, get_f);

  // Psi
  sri::PsiCoreRLE psi_core;
  sdsl::load_from_cache(psi_core, sdsl::conf::KEY_PSI, t_config, true);
  auto psi_select = [&psi_core](auto tt_c, auto tt_rnk) { return psi_core.select(tt_c, tt_rnk); };

  auto get_c = [&alphabet](auto tt_index) { return sri::computeCForSAIndex(alphabet.C, tt_index); };
  auto cumulative = sri::RandomAccessForCRefContainer(std::cref(alphabet.C));

  auto psi = sri::Psi(psi_select, get_c, cumulative);

  // Samples
  sdsl::int_vector<> bwt_run_first; // BWT run heads positions in BWT
  sdsl::load_from_cache(bwt_run_first, conf::KEY_BWT_RUN_FIRST, t_config);

  auto rank_sample = [&bwt_run_first](const auto &tt_k) {
    return std::lower_bound(bwt_run_first.begin(), bwt_run_first.end(), tt_k) - bwt_run_first.begin();
  };

  auto get_link = [&bwt_run_last, &n, &lf, &psi, &rank_sample](const auto &tt_mark_idx) {
    return sri::computeMarkToSampleLinkForPhiForward(bwt_run_last[tt_mark_idx], n, lf, psi, rank_sample);
  };

  // Compute links
  auto [sorted_marks_idx, mark_to_sample_links] = constructMarkToSampleLinks(bwt_run_last_text_pos, get_link);

  sdsl::store_to_cache(sorted_marks_idx, conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_IDX, t_config);
  sdsl::store_to_cache(mark_to_sample_links, conf::KEY_BWT_RUN_LAST_TEXT_POS_SORTED_TO_FIRST_IDX, t_config);
}

inline auto constructMarkToSampleLinksForPhiForwardWithPsiRuns(Config &t_config) {
  using namespace sri::conf;

  sdsl::int_vector<> marks; // Text position of the Psi run tails
  sdsl::load_from_cache(marks, t_config.keys[kPsi][kTail][kTextPos], t_config, true);

  auto get_link = [r = marks.size()](const auto &tt_mark_idx) {
    return (tt_mark_idx + 1) % r;
  };

  auto [sorted_marks_idx, mark_to_sample_links] = constructMarkToSampleLinks(marks, get_link);

  // sdsl::store_to_cache(sorted_marks_idx, t_config.keys[kPsi][kTail][kTextPosAsc][kIdx], t_config);
  sri::store_to_cache(mark_to_sample_links, t_config.keys[kPsi][kTail][kTextPosAsc][kLink], t_config, true);

  return mark_to_sample_links;
}

void constructMarkToSampleLinksForPhiBackward(sdsl::cache_config &t_config) {
  // Marks
  sdsl::int_vector<> marks; // Text position of the first symbol in BWT runs
  sdsl::load_from_cache(marks, conf::KEY_BWT_RUN_FIRST_TEXT_POS, t_config);

  auto get_link = [r = marks.size()](const auto &tt_mark_idx) {
    return (tt_mark_idx + r - 1) % r;
  };

  // Compute links
  auto [sorted_marks_idx, mark_to_sample_links] = constructMarkToSampleLinks(marks, get_link);

  sdsl::store_to_cache(sorted_marks_idx, conf::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_IDX, t_config);
  sdsl::store_to_cache(mark_to_sample_links, conf::KEY_BWT_RUN_FIRST_TEXT_POS_SORTED_TO_LAST_IDX, t_config);
}

template<typename TGetNextMark, typename TGetNextSubmark, typename TReport>
void computeSubmarksValidity(std::size_t t_r_prime,
                             TGetNextMark t_get_next_mark,
                             TGetNextSubmark t_get_next_submark,
                             TReport t_report) {
  // Marks sampled BWT run heads indices in text and if they are trustworthy
  std::size_t mark = t_get_next_mark();
  std::size_t submark = t_get_next_submark();
  for (std::size_t i = 0, j = 0; i < t_r_prime - 1; ++i) {
    std::size_t next_mark = t_get_next_mark();
    std::size_t next_submark = t_get_next_submark();
    if (next_mark != next_submark) {
      // Report current submark as invalid, and what is the next mark to compute valid area
      t_report(i, submark, next_mark);

      do {
        next_mark = t_get_next_mark();
      } while (next_mark != next_submark);
    }

    submark = next_submark;
  }
}

template<typename TValues>
auto constructBitVectorFromIntVector(TValues &t_values, size_t t_bv_size, bool t_init_value) {
  sdsl::bit_vector bv_tmp(t_bv_size, t_init_value);

  for (auto &&item : t_values) {
    bv_tmp[item] = !t_init_value;
  }
  return bv_tmp;
}

template<typename TBitVector, typename TValues, typename TBVRank = typename TBitVector::rank_1_type, typename TBVSelect = typename TBitVector::select_1_type>
void constructBitVectorFromIntVector(TValues &t_values,
                                     const std::string &t_key,
                                     sdsl::cache_config &t_config,
                                     std::size_t t_bv_size,
                                     bool t_init_value) {
  sdsl::bit_vector bv_tmp = constructBitVectorFromIntVector(t_values, t_bv_size, t_init_value);

  TBitVector bv(std::move(bv_tmp));
  sri::store_to_cache(bv, t_key, t_config, true);

  TBVRank bv_rank(&bv);
  sri::store_to_cache(bv_rank, t_key, t_config, true);

  TBVSelect bv_select(&bv);
  sri::store_to_cache(bv_select, t_key, t_config, true);
}

template<typename TBitVector,
  typename TBVRank = typename TBitVector::rank_1_type,
  typename TBVSelect = typename TBitVector::select_1_type>
void constructBitVectorFromIntVector(const std::string& t_key,
                                     sdsl::cache_config& t_config,
                                     std::size_t t_bv_size,
                                     bool t_init_value,
                                     bool t_add_type_hash = false) {
  const auto filename = t_add_type_hash
                          ? sdsl::cache_file_name<sdsl::int_vector<>>(t_key, t_config)
                          : sdsl::cache_file_name(t_key, t_config);
  sdsl::int_vector_buffer<> int_buf(filename);
  constructBitVectorFromIntVector<TBitVector, sdsl::int_vector_buffer<>, TBVRank, TBVSelect>(
    int_buf,
    t_key,
    t_config,
    t_bv_size,
    t_init_value
  );
}

inline void constructSortedIndices(const std::string& t_key,
                                   sdsl::cache_config& t_config,
                                   const std::string& t_out_key,
                                   bool t_add_type_hash = false) {
  sdsl::int_vector<> values;
  sdsl::load_from_cache(values, t_key, t_config, t_add_type_hash);

  auto values_idx = sortIndices(values);

  sdsl::store_to_cache(values_idx, t_out_key, t_config);
}

}

#endif //SRI_CONSTRUCT_H_
