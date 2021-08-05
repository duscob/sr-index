//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/5/21.
//

#ifndef SR_INDEX_PSI_H_
#define SR_INDEX_PSI_H_

namespace sri {

//! Constructs the Psi function using BWT for text.
template<typename TBwt, typename TAlphabet>
auto ConstructPsi(const TBwt &t_bwt, const TAlphabet &t_alphabet) {
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
void ConstructPsi(const TBwt &t_bwt, const TAlphabet &t_alphabet, sdsl::cache_config &t_config) {
  // Store psi
  store_to_cache(ConstructPsi(t_bwt, t_alphabet), sdsl::conf::KEY_PSI, t_config);
}

}

#endif //SR_INDEX_PSI_H_
