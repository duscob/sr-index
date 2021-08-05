//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/27/20.
//

#ifndef RI_R_INDEX_H_
#define RI_R_INDEX_H_

#include <string>
#include <vector>
#include <memory>
#include <experimental/optional>

namespace ri {

class LocateIndex {
 public:
  virtual std::vector<std::size_t> Locate(const std::string &_pattern) const = 0;
};

template<typename TBackwardNav, typename TGetLastValue, typename TComputeAllValues, typename TGetFinalValue>
class RIndex : public LocateIndex {
 public:
  RIndex(const TBackwardNav &t_lf,
         const TGetLastValue &t_get_last_value,
         const TComputeAllValues &t_compute_all_values,
         std::size_t t_bwt_size,
         const TGetFinalValue &t_get_final_sa_value)
      : lf_{t_lf},
        get_last_value_{t_get_last_value},
        compute_all_values_{t_compute_all_values},
        bwt_size_{t_bwt_size},
        get_final_value_{t_get_final_sa_value} {
  }

  std::vector<std::size_t> Locate(const std::string &t_pattern) const override {
    std::vector<std::size_t> values;
    auto report = [&values](const auto &v) { values.emplace_back(v); };

    Locate(t_pattern, report);

    return values;
  }

  template<typename TPattern, typename TReport>
  void Locate(const TPattern &t_pattern, TReport &t_report) const {
    Range range = {0, bwt_size_ - 1};

    auto i = t_pattern.size() - 1;
    auto last_value = get_final_value_(i);

    for (auto it = rbegin(t_pattern); it != rend(t_pattern) && range.first <= range.second; ++it, --i) {
      auto c = *it;

      auto next_range = lf_(range, c);
      last_value = get_last_value_(range, next_range, c, i, last_value);

      range = next_range;
    }

    if (range.first <= range.second) {
      compute_all_values_(range, last_value, t_report);
    }
  }

 private:
  using Range = std::pair<std::size_t, std::size_t>;

  TBackwardNav lf_;
  TGetLastValue get_last_value_;
  TComputeAllValues compute_all_values_;

  std::size_t bwt_size_;
  TGetFinalValue get_final_value_;
};

template<typename TBackwardNav, typename TGetLastValue, typename TComputeAllValues, typename TGetFinalValue>
auto buildRIndex(const TBackwardNav &t_lf,
                 const TGetLastValue &t_get_last_value,
                 const TComputeAllValues &t_compute_all_values,
                 std::size_t t_bwt_size,
                 const TGetFinalValue &t_get_final_sa_value) {
  return RIndex<TBackwardNav, TGetLastValue, TComputeAllValues, TGetFinalValue>(
      t_lf, t_get_last_value, t_compute_all_values, t_bwt_size, t_get_final_sa_value);
}

template<typename TBackwardNav, typename TGetLastValue, typename TComputeAllValues, typename TGetFinalValue>
auto buildSharedPtrRIndex(const TBackwardNav &t_lf,
                          const TGetLastValue &t_get_last_value,
                          const TComputeAllValues &t_compute_all_values,
                          std::size_t t_bwt_size,
                          const TGetFinalValue &t_get_final_sa_value) {
  return std::make_shared<RIndex<TBackwardNav, TGetLastValue, TComputeAllValues, TGetFinalValue>>(
      t_lf, t_get_last_value, t_compute_all_values, t_bwt_size, t_get_final_sa_value);
}

}

#endif //RI_R_INDEX_H_
