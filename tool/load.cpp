#include <sr-index/sr_index.h>

int main(int argc, char *argv[]) {
  std::size_t subsampling_rate = 16;
  sri::SrIndexValidArea<> index(subsampling_rate);

  // Loading the full index from a single file.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Locating occurrences of the given pattern
  std::vector<std::size_t> result = index.Locate("pattern");

  // Processing the result
  // ...

  return 0;
}
