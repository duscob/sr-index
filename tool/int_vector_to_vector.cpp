//
// Created by Dustin Cobas <dustin.cobas@gmail.com> on 2/11/2024.
//

#include <filesystem>
#include <fstream>

#include <gflags/gflags.h>

#include <sdsl/int_vector.hpp>

#include "sr-index/config.h"

DEFINE_string(data, "", "Data file. (MANDATORY)");
DEFINE_string(key, "", "File key. (MANDATORY)");
DEFINE_bool(hash_type, false, "Add hash type to file key.");

int main(int argc, char* argv[]) {
  gflags::AllowCommandLineReparsing();
  gflags::ParseCommandLineFlags(&argc, &argv, false);

  const std::filesystem::path data_path = FLAGS_data;
  sri::Config config(data_path, data_path.parent_path(), sri::SDSL_LIBDIVSUFSORT);

  sdsl::int_vector<> ivec;
  sdsl::load_from_cache(ivec, FLAGS_key, config, FLAGS_hash_type);

  std::ofstream fs(sdsl::cache_file_name(FLAGS_key, config) + ".vec.bin", std::ios::out | std::ios::binary);
  for (sdsl::int_vector<>::value_type item : ivec) {
    fs.write(reinterpret_cast<const char *>(&item), sizeof(sdsl::int_vector<>::value_type));
    // fs.write(&item, sizeof(sdsl::int_vector<>::value_type));
  }
  fs.close();

  return 0;
}
