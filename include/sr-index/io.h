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

namespace sri {

//! Register the existing resource specified by the key to the cache
/*!
 *  \param key        Resource key.
 *  \param config    Cache configuration.
 *
 *  Note: If the resource does not exist under the given key,
 *  it will be not added to the cache configuration.
 */
template <class T>
inline void register_cache_file(std::string const& key, sdsl::cache_config& config) {
  std::string file_name = sdsl::cache_file_name<T>(key, config);
  sdsl::isfstream in(file_name);
  if (in) {  // if file exists, register it.
    config.file_map[key] = file_name;
  }
}

//! Stores the object v as a resource in the cache.
template <class T>
bool store_to_cache(const T& v, const std::string& key, sdsl::cache_config& config, bool add_type_hash = false) {
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

}  // namespace sri

#endif  // SRI_IO_H_
