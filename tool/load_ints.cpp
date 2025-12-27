#include <sr-index/sr_csa.h>

int main(int argc, char* argv[]) {
  std::size_t subsampling_rate = 16;
  sri::SrCSAValidArea<sri::SrCSA<sri::GenericStorage, sri::Alphabet<0>>> index(subsampling_rate);

  // Loading the full index from a single file.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Loading the full index from individual component files
  // If you did not store the full index in a single file,
  // you can use the separated components stored in different files uncommenting the following lines.
  // std::string data_path = "/path/to/data/file";  // only the data filename is required
  // auto output_path = std::filesystem::current_path();
  // uint8_t alphabet_size = 16;  // Input file stores a sequence of 16-bits ints
  // sri::Config config(
  //     data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS, false, alphabet_size, sri::createDefaultKeys<0>());

  // Locating occurrences of the given pattern
  std::vector<uint64_t> pattern = {112, 97, 116, 116, 101, 114, 110};
  std::vector<std::size_t> result = index.Locate(pattern);

  // Processing the result
  // ...

  return 0;
}
