//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/5/21.
//

#ifndef SRI_PSI_H_
#define SRI_PSI_H_

namespace sri {

//! Constructs the Psi function using BWT for text.
template<typename TBwt, typename TAlphabet>
auto constructPsi(TBwt &t_bwt, const TAlphabet &t_alphabet) {
  const auto &n = t_bwt.size();
  const auto &sigma = t_alphabet.sigma;
  sdsl::int_vector<> cnt_chr(sigma, 0, sdsl::bits::hi(n) + 1);
  for (typename TAlphabet::sigma_type i = 0; i < sigma; ++i) {
    cnt_chr[i] = t_alphabet.C[i];
  }

  // Calculate psi
  sdsl::int_vector<> psi(n, 0, sdsl::bits::hi(n) + 1);
  for (std::size_t i = 0; i < n; ++i) {
    psi[cnt_chr[t_alphabet.char2comp[t_bwt[i]]]++] = i;
  }

  return psi;
}

//! Constructs and stores the Psi function using BWT for text.
template<typename TBwt, typename TAlphabet>
void constructPsi(TBwt &t_bwt, const TAlphabet &t_alphabet, sdsl::cache_config &t_config) {
  // Store psi
  store_to_cache(constructPsi(t_bwt, t_alphabet), sdsl::conf::KEY_PSI, t_config);
}

//! Psi function core-representation based on partial psis.
template<typename TBitVector = sdsl::bit_vector, typename TRank = typename TBitVector::rank_1_type, typename TSelect = typename TBitVector::select_1_type>
class PsiCore {
 public:
  const std::vector<TBitVector> &partial_psi = partial_psi_;
  const std::vector<TRank> &rank_partial_psi = rank_partial_psi_;
  const std::vector<TSelect> &select_partial_psi = select_partial_psi_;

  //! Default constructor
  PsiCore() = default;

  //! Constructor
  /**
   * @tparam TCumulativeC Random access container
   * @tparam TPsi Random access container
   * @param t_cumulative_c Cumulative count for the alphabet [0..sigma]
   * @param t_psi Full psi function
   */
  template<typename TCumulativeC, typename TPsi>
  PsiCore(const TCumulativeC &t_cumulative_c, const TPsi &t_psi) {
    auto sigma = t_cumulative_c.size() - 1;

    // Reserve space for psis is critical to avoid invalid pointers in rank/select data structures
    partial_psi_.reserve(sigma);
    rank_partial_psi_.reserve(sigma);
    select_partial_psi_.reserve(sigma);

    auto n = t_cumulative_c[sigma];

    for (std::size_t i = 1; i <= sigma; ++i) {
      // Partial psi for character i
      sdsl::bit_vector psi_c(n, 0);

      // Marks psi values in range of SA corresponding to character i
      for (auto j = t_cumulative_c[i - 1]; j < t_cumulative_c[i]; ++j) {
        psi_c[t_psi[j]] = 1;
      }

      partial_psi_.emplace_back(psi_c);
      rank_partial_psi_.emplace_back(&partial_psi_.back());
      select_partial_psi_.emplace_back(&partial_psi_.back());
    }
  }

  typedef std::size_t size_type;

  //! Serialize method
  size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, const std::string &name = "") const {
    auto child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

    size_type written_bytes = 0;
    written_bytes += sdsl::serialize(partial_psi_, out, child, "m_partial_psi");
    written_bytes += sdsl::serialize(rank_partial_psi_, out, child, "m_rank_partial_psi");
    written_bytes += sdsl::serialize(select_partial_psi_, out, child, "m_select_partial_psi");

    sdsl::structure_tree::add_size(child, written_bytes);

    return written_bytes;
  }

  //! Load method
  void load(std::istream &in) {
    sdsl::load(partial_psi_, in);
    sdsl::load(rank_partial_psi_, in);
    sdsl::load(select_partial_psi_, in);

    for (int i = 0; i < partial_psi_.size(); ++i) {
      rank_partial_psi_[i].set_vector(&partial_psi_[i]);
      select_partial_psi_[i].set_vector(&partial_psi_[i]);
    }
  }

 private:
  // Partial psi function per character. Each bit-vector marks the psi values in the range on SA for the corresponding character.
  std::vector<TBitVector> partial_psi_;
  std::vector<TRank> rank_partial_psi_; // Ranks for partial psis.
  std::vector<TSelect> select_partial_psi_; // Selects for partial psis.
};

//! Gets character c for index in SA using the cumulative counts.
/**
 *
 * @tparam TCumulativeC Random access container
 * @param t_cumulative_c Cumulative count for the alphabet [0..sigma]
 * @param t_index Index in the SA
 * @return Character corresponding to index
 */
template<typename TCumulativeC>
auto computeCForSAIndex(const TCumulativeC &t_cumulative_c, std::size_t t_index) {
  auto upper_bound = std::upper_bound(t_cumulative_c.begin(), t_cumulative_c.end(), t_index);
  return std::distance(t_cumulative_c.begin(), upper_bound) - 1;
}

//! Psi function.
/**
 * @tparam TSelectPartialPsi Select function for partial psis
 * @tparam TGetC Get character function
 * @tparam TCumulativeC Random access container
 */
template<typename TSelectPartialPsi, typename TGetC, typename TCumulativeC>
class Psi {
 public:
  Psi(const TSelectPartialPsi &t_select_partial_psi, const TGetC &t_get_c, const TCumulativeC &t_cumulative_c)
      : select_partial_psi_{t_select_partial_psi}, get_c_{t_get_c}, cumulative_c_(t_cumulative_c) {
  }

  auto operator()(std::size_t t_index) const {
    auto c = get_c_(t_index);
    return select_partial_psi_[c](t_index - cumulative_c_[c] + 1);
  }

 private:
  TSelectPartialPsi select_partial_psi_; // Select function for partial psis
  TGetC get_c_; // Get character function
  TCumulativeC cumulative_c_; // Cumulative count for alphabet [0..sigma]
};

//! LFOnPsi function using partial psi function
/**
 * @tparam TRankPartialPsi Rank function for partial psis
 * @tparam TCumulativeC Random access container
 */
template<typename TRankPartialPsi, typename TCumulativeC>
class LFOnPsi {
 public:
  LFOnPsi(const TRankPartialPsi &t_rank_partial_psi, const TCumulativeC &t_cumulative_c)
      : rank_partial_psi_{t_rank_partial_psi}, cumulative_c_{t_cumulative_c} {
  }

  //! LFOnPsi function
  /**
   * @tparam TRange Range type {sp; ep}
   * @tparam TChar Character type for compact alphabet
   * @param t_range Range [sp; ep]
   * @param t_c Character for compact alphabet in [0..sigma]
   * @return [new_sp; new_ep]
   */
  template<typename TRange, typename TChar>
  auto operator()(const TRange &t_range, const TChar &t_c) const {
    auto[sp, ep] = t_range;

    // Number of c before the interval
    auto c_before_sp = rank_partial_psi_[t_c](sp);

    // Number of c before the interval + number of c inside the interval range
    auto c_until_ep = rank_partial_psi_[t_c](ep + 1);

    // If there are no c in the interval, return empty range
    if (c_before_sp == c_until_ep)
      return TRange{1, 0};

    // Number of characters smaller than c
    auto prev_to_c = cumulative_c_[t_c];

    return TRange{prev_to_c + c_before_sp, prev_to_c + c_until_ep - 1};
  }

 private:
  TRankPartialPsi rank_partial_psi_; // Rank function for partial psis
  TCumulativeC cumulative_c_; // Cumulative count for alphabet [0..sigma]
};
}

#endif //SRI_PSI_H_
