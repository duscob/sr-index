//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/24/21.
//

#ifndef SRI_IO_H_
#define SRI_IO_H_

#include <iostream>
#include <string>
#include <utility>

#include <sdsl/config.hpp>
#include <sdsl/io.hpp>
#include <sdsl/util.hpp>

namespace sri {

template <typename T, typename = void>
struct is_serializable : std::false_type {};

template <typename T>
struct is_serializable<T,
                       std::void_t<decltype(serialize(std::declval<T>(),
                                                      std::declval<std::ostream&>(),
                                                      std::declval<sdsl::structure_tree_node*>(),
                                                      std::declval<const std::string&>()))>> : std::true_type {};

template <typename T, typename = void>
struct is_loadable : std::false_type {};

template <typename T>
struct is_loadable<T, std::void_t<decltype(load(std::declval<T>(), std::declval<std::istream&>()))>> : std::true_type {
};

}  // namespace sri

namespace std {

template <typename X, typename Y>
std::enable_if_t<!sri::is_serializable<std::pair<X, Y>>::value, uint64_t> serialize(
    const std::pair<X, Y>& x,
    std::ostream& out,
    sdsl::structure_tree_node* v = nullptr,
    const std::string& name = "") {
  return serialize(x.first, out, v, name) + serialize(x.second, out, v, name);
}

template <typename X, typename Y>
std::enable_if_t<!sri::is_loadable<std::pair<X, Y>>::value, void> load(std::pair<X, Y>& x, std::istream& in) {
  using sdsl::load;
  load(x.first, in);
  load(x.second, in);
}

}  // namespace std


#endif  // SRI_IO_H_
