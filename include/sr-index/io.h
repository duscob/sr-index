//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/24/21.
//

#ifndef SRI_IO_H_
#define SRI_IO_H_

#include <utility>
#include <iostream>
#include <string>

#include "sdsl/io.hpp"

namespace std {

template<typename X, typename Y>
uint64_t serialize(const std::pair<X, Y> &x,
                   std::ostream &out,
                   sdsl::structure_tree_node *v = nullptr,
                   const std::string &name = "") {
  return serialize(x.first, out, v, name) + serialize(x.second, out, v, name);
}

template<typename X, typename Y>
void load(std::pair<X, Y> &x, std::istream &in) {
  using sdsl::load;
  load(x.first, in);
  load(x.second, in);
}

}

#endif //SRI_IO_H_
