#include <sr-index/r_csa.h>

int main(int argc, char* argv[]) {
  std::string data_path = "/path/to/data/file";
  auto output_path = std::filesystem::current_path();
  uint8_t alphabet_size = 16;  // Input file stores a sequence of 16-bits ints
  sri::Config config(
      data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS, false, alphabet_size, sri::createDefaultKeys<0>());

  sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> index;

  // Constructs only the count-required components (alphabet + psi-RLE), skipping the locate
  // artifacts. After this call the index supports `Count` queries.
  sri::construct(index, data_path, config);

  // The serialized file is byte-identical to the prefix of an `RCSA` locate bundle.
  std::string index_file = "/path/to/index/file";
  sdsl::store_to_file(index, index_file);

  return 0;
}
