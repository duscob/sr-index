//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/24/21.
//

#ifndef SRI_CODER_H_
#define SRI_CODER_H_

#include <sdsl/coder.hpp>
#include <sdsl/bits.hpp>

namespace sri {

template<typename T>
auto decode(const uint64_t *&data, std::size_t &start_idx);

template<>
auto decode<sdsl::coder::elias_delta>(const uint64_t *&d, std::size_t &start_idx) {
  d += (start_idx >> 6);
  uint64_t value = 0;
  uint8_t offset = start_idx & 0x3F;

  auto len_1_len = sdsl::bits::read_unary_and_move(d, offset); // read length of length of x
  if (!len_1_len) {
    value += 1;
  } else {
    auto len = sdsl::bits::read_int_and_move(d, offset, len_1_len) + (1ULL << len_1_len);
    value += sdsl::bits::read_int_and_move(d, offset, len - 1) + (len - 1 < 64) * (1ULL << (len - 1));
  }

  start_idx = offset;
  return value;
}

}

#endif //SRI_CODER_H_
