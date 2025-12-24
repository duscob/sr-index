//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 2/12/22.
//

#ifndef SRI_ALPHABET_H_
#define SRI_ALPHABET_H_

#include <cstdint>

#include <sdsl/csa_alphabet_strategy.hpp>

namespace sri {

template <uint8_t t_width>
struct alphabet_trait {
  typedef sdsl::int_alphabet<> type;
};

template <>
struct alphabet_trait<8> {
  typedef sdsl::byte_alphabet type;
};

template <uint8_t t_width = 8>
class Alphabet : public alphabet_trait<t_width>::type {};

}  // namespace sri

#endif  // SRI_ALPHABET_H_
