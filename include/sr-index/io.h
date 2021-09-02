//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 8/24/21.
//

#ifndef SRI_IO_H_
#define SRI_IO_H_

#include <utility>
#include <iostream>
#include <string>

#include <sdsl/io.hpp>
#include <sdsl/config.hpp>
#include <sdsl/util.hpp>

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

namespace sri {

//! Stores the object v as a resource in the cache.
template<class T>
bool store_to_cache(const T &v, const std::string &key, sdsl::cache_config &config, bool add_type_hash = false) {
  std::string file;
  if (add_type_hash) {
    file = sdsl::cache_file_name<T>(key, config);
  } else {
    file = sdsl::cache_file_name(key, config);
  }
  if (sdsl::store_to_file(v, file)) {
    config.file_map[key + (add_type_hash ? "_" + sdsl::util::class_to_hash(T()) : "")] = file;
    return true;
  } else {
    std::cerr << "WARNING: store_to_cache: could not store file `" << file << "`" << std::endl;
    return false;
  }
}

}

#endif //SRI_IO_H_
