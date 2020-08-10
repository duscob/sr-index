//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/8/20.
//

#ifndef RI_BWT_H_
#define RI_BWT_H_

#include <cstddef>
#include <cassert>

namespace ri {

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
  for (int i = 1; i < last_pos; ++i) {
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

  // Last BWT value
  text_pos = get_text_pos(last_pos);
  bwt_symbol = t_get_text_at(text_pos);
  t_report_bwt(last_pos, bwt_symbol);

  // Last position ends the last BWT run
  t_report_bwt_run_last(n_runs, last_pos, bwt_symbol, text_pos);

  return n_runs + 1;
}

}

#endif //RI_BWT_H_
