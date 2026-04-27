#include <sr-index/r_index.h>

int main(int argc, char* argv[]) {
  sri::RIndexCount<sri::GenericStorage, sri::Alphabet<8>> index;

  // Loading from a single file: works whether the file was built by `RIndexCount` (count-only
  // bundle) or by `RIndex` (locate bundle). In the latter case the count-only loader reads only
  // the prefix and stops.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Alternatively, load from individual cached component files:
  // std::string data_path = "/path/to/data/file";
  // auto output_path = std::filesystem::current_path();
  // sri::Config config(data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS);
  // index.load(config);

  // Counting occurrences of the given pattern (returns the BWT range [sp, ep]).
  auto [sp, ep] = index.Count("pattern");
  std::size_t occurrences = (sp <= ep) ? (ep - sp + 1) : 0;
  (void)occurrences;

  return 0;
}
