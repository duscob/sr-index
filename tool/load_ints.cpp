#include <sr-index/sr_index.h>
#include <sr_csa.h>

int main(int argc, char* argv[]) {
  std::size_t subsampling_rate = 16;
  sri::SrIndexValidArea<sri::GenericStorage, sri::Alphabet<0>> index(subsampling_rate);

  // Loading the full index from a single file.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Locating occurrences of the given pattern
  std::vector<uint64_t> pattern = {112, 97, 116, 116, 101, 114, 110};
  std::vector<std::size_t> result = index.Locate(pattern);

  // Processing the result
  // ...

  return 0;
}
