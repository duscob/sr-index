//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/19/2023.
//

#ifndef SRI_CONFIG_H_
#define SRI_CONFIG_H_

#include <filesystem>

#include <sdsl/config.hpp>

namespace sri {

struct Config : public sdsl::cache_config {
  std::filesystem::path data_path;

  Config() = default;

  Config(const std::filesystem::path &t_data_path, const std::filesystem::path &t_output_dir)
      : data_path(t_data_path),
        cache_config(false, t_output_dir, t_data_path.filename()) {
  }
};

} // namespace sri

#endif //SRI_CONFIG_H_
