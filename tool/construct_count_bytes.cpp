#include <sr-index/r_index.h>

int main(int argc, char* argv[]) {
  std::string data_path = "/path/to/data/file";
  auto output_path = std::filesystem::current_path();
  sri::Config config(data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS);

  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> index;

  // Constructs only the count-required components (alphabet + bwt-rle), skipping the locate
  // artifacts. After this call the index supports `Count` queries.
  sri::construct(index, data_path, config);

  // Serializing the count-only index. The resulting file is byte-identical to the prefix of a
  // locate bundle built with `RIndex` over the same data — a count-only loader can read either.
  std::string index_file = "/path/to/index/file";
  sdsl::store_to_file(index, index_file);

  return 0;
}
