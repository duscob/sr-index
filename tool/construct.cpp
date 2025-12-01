#include <sr-index/sr_index.h>

int main(int argc, char *argv[]) {
  std::string data_path = "/path/to/data/file";
  auto output_path = std::filesystem::current_path();
  sri::Config config(data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS);

  std::size_t subsampling_rate = 16;
  sri::SrIndexValidArea<> index(subsampling_rate);

  // Constructing required components and serializing them in separated files.
  // After construction, the full index is loaded in the `index` variable
  // and is ready to be used with locate and count operations.
  sri::construct(index, data_path, config);

  // Serializing the full index in a single file.
  // This is not required if you keep all the files created in the output_path.
  std::string index_file = "/path/to/index/file";
  sdsl::store_to_file(index, index_file);

  return 0;
}
