//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/19/2023.
//

#ifndef SRI_CONFIG_H_
#define SRI_CONFIG_H_

#include <filesystem>
#include <utility>

#include <sdsl/config.hpp>

namespace sri {

struct Config : public sdsl::cache_config {
  std::filesystem::path data_path;

  explicit Config(std::filesystem::path t_data_path = "./data")
      : data_path(std::move(t_data_path)),
        cache_config(false, t_data_path.parent_path(), t_data_path.filename()) {
  }

  Config(const std::filesystem::path &t_data_dir, const std::filesystem::path &t_data_name)
      : data_path(t_data_dir / t_data_name),
        cache_config(false, t_data_dir, t_data_name) {
  }
};

} // namespace sri

#endif //SRI_CONFIG_H_
