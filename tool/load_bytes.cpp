#include <sr-index/sr_index.h>

int main(int argc, char* argv[]) {
  std::size_t subsampling_rate = 16;
  sri::SrIndexValidArea<sri::GenericStorage, sri::Alphabet<8>> index(subsampling_rate);

  // Loading the full index from a single file.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Loading the full index from individual component files
  // If you did not store the full index in a single file,
  // you can use the separated components stored in different files uncommenting the following lines.
  // std::string data_path = "/path/to/data/file";  // only the data filename is required
  // auto output_path = std::filesystem::current_path();
  // sri::Config config(data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS);
  // index.load(config);


  // Locating occurrences of the given pattern
  std::vector<std::size_t> result = index.Locate("pattern");

  // Processing the result
  // ...

  return 0;
}
