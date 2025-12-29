$sr$-index and $sr$-csa: Fast and Small Indexes for Repetitive Texts
=====


> Dustin Cobas, Travis Gagie, and Gonzalo Navarro. Fast and Small Subsampled R-indexes. ACM Transactions on Algorithms,
> Volume 22, Issue 1, Article 7 (October 2025), Pages 1–39.

Overview
-----

This repository provides a C++ implementation of our **$sr$-index**, a variant of $r$-index that limits a large fraction
of the space to $\mathcal{O}(\min(r, n/s))$ for a text of length $n$ with $r$ runs in its Burrows–Wheeler Transform (
BWT) and a given parameter $s$, at the expense of multiplying by $s$ the time per occurrence reported.
The **$sr$-index** is obtained by carefully subsampling the text positions indexed by the $r$-index, in a way that we
prove is still able to support pattern matching with guaranteed performance.

Implementations of our **$r$-csa** and **$sr$-csa** indexes are provided too. Just like the $r$-index adapts the
well-known FM-Index to repetitive texts, the **$r$-csa** adapts Sadakane’s Compressed Suffix Array (CSA) to this case.
The **$sr$-csa** is the corresponding subsampled version of the **$r$-csa**.

Our experiments demonstrate that the theoretical analysis falls short in describing the practical advantages of the *
*$sr$-index** and **$sr$-csa**. For example, the **$sr$-index** retains the performance of the $r$-index while using
1.5–4.0 times less space, sharply outperforming virtually every other compressed index on repetitive texts in both time
and space.

***


Implemented Indexes
-----

This repository provides implementations of the following indexes:

- **$r$-index** (`RIndex`): Base $r$-index for repetitive texts, using the BWT runs structure and $\phi$ function.
- **$sr$-index** (`SrIndex`): Subsampled version of the $r$-index with reduced space usage.
- **$sr$-index with valid marks** (`SrIndexValidMark`): Enhanced $sr$-index that marks valid remaining samplings
  for $\phi$ function.
- **$sr$-index with valid areas** (`SrIndexValidArea`): Enhanced $sr$-index that tracks valid text areas of remaing
  samplings for $\phi$ function.
- **$r$-csa** (`RCSA`): Compressed Suffix Array adapted for repetitive texts, analogous to the $r$-index (but
  using $\psi$ array instead of BWT, and $\phi^{-1}$ function)
- **$sr$-csa** (`SrCSA`): Subsampled version of the $r$-csa with reduced space usage
- **$sr$-csa with valid marks** (`SrCSAValidMark`): Enhanced $sr$-csa that marks valid remaining samplings
  for $\phi^{-1}$ function.
- **$sr$-csa with valid areas** (`SrCSAValidArea`): Enhanced $sr$-csa that tracks valid text areas of remaing samplings
  for $\phi^{-1}$ function.

All indexes support:

- **Count**: Count the number of occurrences of a pattern.
- **Locate**: Find all occurrences of a pattern in the indexed text.

The indexes work with both byte alphabets (text files) and integer alphabets (binary files) through generic template
parameters.

***


Requirements
-----

- C++17-capable compiler
    - GCC 8+ (note: special linking for std::filesystem on GCC < 9 is handled in CMake)
    - Clang 7+ (or compatible)
- CMake 3.10+ (Build system and dependency management)
- Git and internet access during the first configure to fetch dependencies
- POSIX-like environment recommended; Windows may work but is not a primary target

### Third‑party Dependencies

All core dependencies are fetched and built automatically by CMake:

- [SDSL (Succinct Data Structure Library)](https://github.com/duscob/sdsl-lite.git)
  #fork [duscob](https://github.com/duscob)
- [nlohmann/json](https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz) #v3.12.0 (tarball)
- [Big-BWT](https://github.com/duscob/Big-BWT.git) #fork [duscob](https://github.com/duscob): Optional
- [gflags](https://github.com/gflags/gflags.git) #v2.3.0: Optional (CLI flags for tools and benchmarks)
- [Google Test](https://github.com/google/googletest.git) #release-1.11.0: Optional (for tests)
- [Google Benchmark](https://github.com/google/benchmark.git) #v1.9.4: Optional (for benchmarks)

***


Library Consumption: Adding **$sr$-indexes** into an Existing CMake Project
-----

This repository provides a header-only C++ library that implements **$sr$-index** and **$sr$-csa**.


#### Fetching and Adding **$sr$-indexes**

The following snippet adds the **$sr$-index** library to an existing CMake project using FetchContent.

```cmake
set(ExternalProjectName sr-index)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/sr-index.git
        GIT_TAG feature/cmake
        FIND_PACKAGE_ARGS
)

set(SR-INDEX_ENABLE_TESTS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_TOOLS OFF CACHE BOOL "")
set(SR-INDEX_ENABLE_BENCHMARKS OFF CACHE BOOL "")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)
```

### Byte Alphabet (Text File)

#### Constructing and Storing **$sr$-index**

The following snippet constructs a $sr$-index (with valid areas) for a given text file and stores it in a given
directory.

```c++
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
```

#### Loading and Using **$sr$-index**

The following snippet loads a $sr$-index (with valid areas) from a single file and uses it to locate occurrences of a
given pattern.

```c++
#include <sr-index/sr_index.h>


int main(int argc, char *argv[]) {
  std::size_t subsampling_rate = 16;
  sri::SrIndexValidArea<> index(subsampling_rate);

  // Loading the full index from a single file.
  std::string index_file = "/path/to/index/file";
  sdsl::load_from_file(index, index_file);
    
  // Loading the full index from individual component files
  // If you did not store the full index in a single file,
  // you can use the separated components stores in different files uncommenting the following lines.
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

```

### Int Alphabet (Binary File)

#### Constructing and Storing **$sr$-csa**

The following snippet constructs a $sr$-csa (with valid areas) for a given binary file of 16-bit ints and stores it in a
given directory.

```c++
#include <sr-index/sr_csa.h>

int main(int argc, char* argv[]) {
  std::string data_path = "/path/to/data/file";
  auto output_path = std::filesystem::current_path();
  uint8_t alphabet_size = 16;  // Input file stores a sequence of 16-bits ints
  sri::Config config(
      data_path, output_path, sri::SAAlgo::SDSL_SE_SAIS, false, alphabet_size, sri::createDefaultKeys<0>());

  std::size_t subsampling_rate = 16;
  sri::SrCSAValidArea<sri::SrCSA<sri::GenericStorage, sri::Alphabet<0>>> index(subsampling_rate);

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
```

#### Loading and Using **$sr$-csa**

The following snippet loads a $sr$-csa (with valid areas) from a single file and uses it to locate occurrences of a
given pattern.

```c++
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
```

***


Standalone Build and Use
-----

Basic out-of-source build (Release):

```bash
git clone https://github.com/duscob/sr-index.git
cd sr-index
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DSR-INDEX_ENABLE_TESTS=ON \
  -DSR-INDEX_ENABLE_TOOLS=ON \
  -DSR-INDEX_ENABLE_BENCHMARKS=ON
cmake --build build --config Release --parallel
```

Configuration options:

- `SR-INDEX_ENABLE_TESTS` (ON/OFF, default ON): build unit tests
- `SR-INDEX_ENABLE_TOOLS` (ON/OFF, default ON): build auxiliary tools
- `SR-INDEX_ENABLE_BENCHMARKS` (ON/OFF, default ON): build benchmarks
- `SR-INDEX_ENABLE_BIG_BWT` (ON/OFF, default OFF): enable Big-BWT support and define `SRI_USE_BIG_BWT=1`

This project builds an interface library target named `sr-index` and optional executables. There is no install step
provided by default. If you need installation, consider adding `install()` rules or consuming the project via
`add_subdirectory`.

### Entry points (executables)

#### Tools (enabled with `SR-INDEX_ENABLE_TOOLS=ON`):

- `int_vector_to_vector` — converts an `sdsl::int_vector<>` cached object into a raw binary vector file
    - Flags (gflags):
        - `--data` (string, required): path to the data file (used to derive cache config)
        - `--key` (string, required): cache key
        - `--hash_type` (bool): whether to incorporate hash type into the cache key
    - Example:
      ```bash
      ./build/int_vector_to_vector --data=corpus.txt --key=psi_core_bv --hash_type=false
      ```

#### Benchmarks (enabled with `SR-INDEX_ENABLE_BENCHMARKS=ON`):

- Top-level:
    - `bm_build_ri_items`
    - `bm_locate`
- **$sr$-index** benchmarks (`benchmark/sr-index`):
    - `bm_construct_ri`, `bm_locate_ri`, `bm_count_ri`
- **$sr$-csa** benchmarks (`benchmark/sr-csa`):
    - `bm_construct_csa`, `bm_locate_csa`, `bm_count_csa`

These are Google Benchmark executables; some may also accept gflags. Usage varies by benchmark.

TODO: Document specific CLI flags and input formats for each benchmark once stabilized.


### Running tests

Build with tests enabled, then use CTest:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSR-INDEX_ENABLE_TESTS=ON
cmake --build build --parallel
cd build
ctest --output-on-failure
```

Environment variables and external tools
-----

### Big-BWT

Big-BWT is an optional dependency that can be enabled at build time using the `SR-INDEX_ENABLE_BIG_BWT` CMake option:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSR-INDEX_ENABLE_BIG_BWT=ON
cmake --build build --parallel
```

When enabled:
- The `SRI_USE_BIG_BWT=1` preprocessor definition is set, enabling Big-BWT code paths
- During index construction, the code may invoke an external Big-BWT executable to compute the BWT
- The Big-BWT executable path is configured at build time via the `-DBIGBWT_EXE` definition, which CMake sets to the fetched Big-BWT binary directory
- No runtime environment variables are currently required

If you need to use a system-provided Big-BWT instead of the fetched one, you may adjust CMake to point `BIGBWT_EXE`
accordingly. TODO: Provide an option to override `BIGBWT_EXE` via CMake.

***


Project structure
-----

```
include/                 # Header-only library
  sr-index/
    *.h, *.hpp           # Core data structures and algorithms
benchmark/               # Google Benchmark sources (built when benchmarks are enabled)
  sr-index/*             # r-index and sr-index benchmarks
  sr-csa/*               # r-csa and sr-csa benchmarks
tool/                    # Small standalone tools (built when tools are enabled)
  int_vector_to_vector.cpp
test/                    # GoogleTest unit tests (built when tests are enabled)
  *.cpp                  # GoogleTest unit tests
cmake/                   # CMake helper modules and dependency configs
  *.cmake                # Dependency and helper modules
CMakeLists.txt           # Root CMake configuration
```

Key headers (non-exhaustive): `r_index.h`, `sr_index.h`, `r_csa.h`, `sr_csa.h`, `lf.h`, `bwt.h`, `psi.h`, `phi.h`,
`toehold.h`, `sampling.h`, `rle_string.hpp`, `sparse_sd_vector.hpp`, `sparse_hyb_vector.hpp`, `huff_string.hpp`,
`construct*.h(pp)`.

### Development notes

- Compiler flags are configured for Debug/Release/RelWithDebInfo in `CMakeLists.txt`.
- On GCC 8, `-lstdc++fs` is linked automatically for std::filesystem support.

***


License
-----

TODO: Add license file and specify the licensing terms for this repository.

***


Citation
-----

If you use this code in your research, please cite our paper:

```bibtex
@article{10.1145/3750729,
    author = {Cobas, Dustin and Gagie, Travis and Navarro, Gonzalo},
    title = {Fast and Small Subsampled R-indexes},
    year = {2025},
    issue_date = {January 2026},
    volume = {22},
    number = {1},
    issn = {1549-6325},
    url = {https://doi.org/10.1145/3750729},
    doi = {10.1145/3750729},
    journal = {ACM Trans. Algorithms},
    month = oct,
    articleno = {7},
    numpages = {29},
}
```
