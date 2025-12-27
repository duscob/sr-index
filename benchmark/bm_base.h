//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 12/25/25.
//

#pragma once

#include <functional>

#include "sr-index/alphabet.h"
#include "sr-index/index_base.h"


using ExternalGenericStorage = std::reference_wrapper<sri::GenericStorage>;

using Sequence = sri::Alphabet<>::string_type;

template <typename Char>
std::ostream& operator<<(std::ostream& os, const std::vector<Char>& seq) {
  for (auto c : seq) {
    os << c << " ";
  }

  return os;
}

template <typename Char>
std::vector<Char> base64_decode(const std::vector<Char>& seq, bool remove_linebreaks = false) {
  return seq;
}
