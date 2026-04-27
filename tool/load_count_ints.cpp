#include <sr-index/r_csa.h>

int main(int argc, char* argv[]) {
  sri::RCSACount<sri::GenericStorage, sri::Alphabet<0>> index;

  // Loading from a single file: works whether the file was built by `RCSACount` (count-only
  // bundle) or by `RCSA` (locate bundle). In the latter case the count-only loader reads only
  // the prefix and stops.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);

  // Counting occurrences of the given pattern (returns the BWT range [sp, ep]).
  std::vector<uint64_t> pattern = {112, 97, 116, 116, 101, 114, 110};
  auto [sp, ep] = index.Count(pattern);
  std::size_t occurrences = (sp <= ep) ? (ep - sp + 1) : 0;
  (void)occurrences;

  return 0;
}
