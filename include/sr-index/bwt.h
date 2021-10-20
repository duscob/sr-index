//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/8/20.
//

#ifndef RI_BWT_H_
#define RI_BWT_H_

#include <cstddef>
#include <cassert>
#include <memory>

namespace sri {

template<typename TGetSAAt, typename TGetTextAt, typename TReportBWT, typename TReportBWTRunFirst, typename TReportBWTRunLast>
std::size_t computeBWT(std::size_t t_n,
                       const TGetSAAt &t_get_sa_at,
                       const TGetTextAt &t_get_text_at,
                       const TReportBWT &t_report_bwt,
                       const TReportBWTRunFirst &t_report_bwt_run_first,
                       const TReportBWTRunLast &t_report_bwt_run_last) {
  assert(0 < t_n);

  const auto last_pos = t_n - 1;

  auto get_text_pos = [&t_get_sa_at, &last_pos](auto bwt_idx) {
    auto next_pos = t_get_sa_at(bwt_idx);
    return 0 < next_pos ? next_pos - 1 : last_pos;
  };

  size_t n_runs = 0; // # BWT runs

  // First BWT value
  auto text_pos = get_text_pos(0);
  auto bwt_symbol = t_get_text_at(text_pos);
  t_report_bwt(0, bwt_symbol);

  // First position starts the first BWT run.
  t_report_bwt_run_first(n_runs, 0, bwt_symbol, text_pos);

  auto prev_bwt_symbol = bwt_symbol;
  auto prev_text_pos = text_pos;
  for (int i = 1; i <= last_pos; ++i) {
    text_pos = get_text_pos(i);
    bwt_symbol = t_get_text_at(text_pos);
    t_report_bwt(i, bwt_symbol);

    if (bwt_symbol != prev_bwt_symbol) {
      // Last position of the previous BWT run
      t_report_bwt_run_last(n_runs, i - 1, prev_bwt_symbol, prev_text_pos);

      ++n_runs;

      // First position of the current BWT run
      t_report_bwt_run_first(n_runs, i, bwt_symbol, text_pos);

      prev_bwt_symbol = bwt_symbol;
    }

    prev_text_pos = text_pos;
  }

  // Last position ends the last BWT run
  t_report_bwt_run_last(n_runs, last_pos, bwt_symbol, text_pos);

  return n_runs + 1;
}

/// Backward navigation (Last to First) on BWT
template<typename TGetCRankOnBWT, typename TGetF, typename Range, typename TChar>
Range computeLF(const TGetCRankOnBWT &get_c_rank_on_bwt,
                const TGetF &get_f,
                const Range &range,
                TChar c,
                std::size_t bwt_size,
                TChar max_c = 255) {
  auto prev_to_c = get_f(c);
  //if character does not appear in the text, return empty pair
  if ((c == max_c and prev_to_c == bwt_size) || prev_to_c >= get_f(c + 1))
    return {1, 0};

  //number of c before the interval
  auto c_before = get_c_rank_on_bwt(range.first, c);

  //number of c inside the interval range
  auto c_inside = get_c_rank_on_bwt(range.second + 1, c) - c_before;

  //if there are no c in the interval, return empty range
  if (c_inside == 0)
    return {1, 0};

  auto l = prev_to_c + c_before;

  return {l, l + c_inside - 1};
}

/// Backward navigation on BWT
template<typename TGetCRankOnBWT, typename TGetF, typename TChar = unsigned char>
class LFOnBWT {
 public:
  using Range = std::pair<std::size_t, std::size_t>;

  LFOnBWT(const TGetCRankOnBWT &t_get_c_rank_on_bwt,
          const TGetF &t_f,
          std::size_t t_bwt_size,
          TChar t_max_c = 255)
      : get_c_rank_on_bwt_{t_get_c_rank_on_bwt}, f_{t_f}, bwt_size_{t_bwt_size}, max_c_{t_max_c} {
  }

  Range operator()(const Range &range, TChar c) const {
    return computeLF(get_c_rank_on_bwt_, f_, range, c, bwt_size_, max_c_);
  }

 private:
  TGetCRankOnBWT get_c_rank_on_bwt_; // BWT rank for a given char and position
  TGetF f_; // Accumulative frequencies by symbol

  std::size_t bwt_size_; // Size of the BWT
  TChar max_c_; // Maximum symbol
};

template<typename TGetCRankOnBWT, typename TGetF, typename TChar = unsigned char>
auto buildLF(const TGetCRankOnBWT &t_get_c_rank_on_bwt,
             const TGetF &t_f,
             std::size_t t_bwt_size,
             TChar t_max_c = 255) {
  return LFOnBWT<TGetCRankOnBWT, TGetF, TChar>(t_get_c_rank_on_bwt, t_f, t_bwt_size, t_max_c);
}

template<typename TGetChar, typename TGetCRankOnBWT, typename TGetF>
class BasicLF {
 public:
  BasicLF(const TGetChar &t_get_char, const TGetCRankOnBWT &t_get_c_rank_on_bwt, const TGetF &t_f)
      : get_char_{t_get_char}, get_c_rank_on_bwt_{t_get_c_rank_on_bwt}, f_(t_f) {
  }

  std::size_t operator()(std::size_t t_pos) const {
    auto c = get_char_(t_pos);
    return f_(c) + get_c_rank_on_bwt_(t_pos, c);
  }

 private:
  TGetChar get_char_;
  TGetCRankOnBWT get_c_rank_on_bwt_;

  TGetF f_;
};

template<typename TGetChar, typename TGetCRankOnBWT, typename TGetF>
auto buildBasicLF(const TGetChar &t_get_char, const TGetCRankOnBWT &t_get_c_rank_on_bwt, const TGetF &t_f) {
  return BasicLF<TGetChar, TGetCRankOnBWT, TGetF>(t_get_char, t_get_c_rank_on_bwt, t_f);
}

}

#endif //RI_BWT_H_
