//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_R_INDEX_H_
#define RI_R_INDEX_H_

#include <string>
#include <vector>
#include <memory>

namespace sri {

class LocateIndex {
 public:
  virtual std::vector<std::size_t> Locate(const std::string &_pattern) const = 0;
};

template<typename TBackwardNav, typename TGetLastValue, typename TComputeAllValues, typename TGetFinalValue, typename TGetSymbol, typename TCreateFullRange, typename TIsRangeEmpty>
class RIndex : public LocateIndex {
 public:
  RIndex(const TBackwardNav &t_lf,
         const TGetLastValue &t_get_last_value,
         const TComputeAllValues &t_compute_all_values,
         std::size_t t_bwt_size,
         const TGetFinalValue &t_get_final_sa_value,
         const TGetSymbol &t_get_symbol,
         const TCreateFullRange &t_create_full_range,
         const TIsRangeEmpty &t_is_range_empty)
      : lf_{t_lf},
        get_last_value_{t_get_last_value},
        compute_all_values_{t_compute_all_values},
        bwt_size_{t_bwt_size},
        get_final_value_{t_get_final_sa_value},
        get_symbol_{t_get_symbol},
        create_full_range_{t_create_full_range},
        is_range_empty_{t_is_range_empty} {
  }

  std::vector<std::size_t> Locate(const std::string &t_pattern) const override {
    std::vector<std::size_t> values;
    auto report = [&values](const auto &v) { values.emplace_back(v); };

    Locate(t_pattern, report);

    return values;
  }

  template<typename TPattern, typename TReport>
  void Locate(const TPattern &t_pattern, TReport &t_report) const {
    auto range = create_full_range_(bwt_size_);

    auto i = t_pattern.size() - 1;
    auto last_value = get_final_value_(i); //TODO use default value (step == 0) instead of get_final_value_

    for (auto it = rbegin(t_pattern); it != rend(t_pattern) && !is_range_empty_(range); ++it, --i) {
      auto c = get_symbol_(*it);

      auto next_range = lf_(range, c);
      // TODO refactor get_last_value to avoid unneeded copy of last_value in trivial case (parameter as reference)
      last_value = get_last_value_(range, next_range, c, i, last_value);

      range = next_range;
    }

    if (!is_range_empty_(range)) {
      compute_all_values_(range, last_value, t_report);
    }
  }

 private:

  TBackwardNav lf_;
  TGetLastValue get_last_value_;
  TComputeAllValues compute_all_values_;

  std::size_t bwt_size_;
  TGetFinalValue get_final_value_;

  TGetSymbol get_symbol_;

  TCreateFullRange create_full_range_;
  TIsRangeEmpty is_range_empty_;
};


template<typename TBackwardNav, typename TGetLastValue, typename TComputeAllValues, typename TGetFinalValue, typename TGetSymbol>
auto buildSharedPtrRIndex(const TBackwardNav &t_lf,
                          const TGetLastValue &t_get_last_value,
                          const TComputeAllValues &t_compute_all_values,
                          std::size_t t_bwt_size,
                          const TGetFinalValue &t_get_final_sa_value,
                          const TGetSymbol &t_get_symbol) {
  using Range = std::pair<std::size_t, std::size_t>;
  using TFnCreateFullRange = std::function<Range(std::size_t)>;
  auto create_full_range = [](auto tt_seq_size) { return Range{0, tt_seq_size - 1}; };

  using TFnIsRangeEmpty = std::function<bool(const Range&)>;
  auto is_range_empty = [](const auto &tt_range) { return  tt_range.second < tt_range.first; };

  return std::make_shared<RIndex<TBackwardNav, TGetLastValue, TComputeAllValues, TGetFinalValue, TGetSymbol, TFnCreateFullRange, TFnIsRangeEmpty>>(
      t_lf, t_get_last_value, t_compute_all_values, t_bwt_size, t_get_final_sa_value, t_get_symbol, create_full_range, is_range_empty);
}

}

#endif //RI_R_INDEX_H_
