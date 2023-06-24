//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 6/19/2023.
//

#ifndef SRI_CONFIG_H_
#define SRI_CONFIG_H_

#include <filesystem>

#include <sdsl/config.hpp>
#include <sdsl/construct_config.hpp>

namespace sri {

enum SAAlgo {
  SDSL_LIBDIVSUFSORT,
  SDSL_SE_SAIS,
  BIG_BWT
};

struct Config : public sdsl::cache_config {
  std::filesystem::path data_path;
  SAAlgo sa_algo = SDSL_LIBDIVSUFSORT;

  Config() = default;

  Config(const std::filesystem::path &t_data_path,
         const std::filesystem::path &t_output_dir,
         SAAlgo t_sa_algo = SDSL_LIBDIVSUFSORT)
      : data_path(t_data_path),
        sa_algo(t_sa_algo),
        cache_config(false, t_output_dir, t_data_path.filename()) {
  }
};

} // namespace sri

#endif //SRI_CONFIG_H_
