option(sr-index_build_tests "Build all of sr-index's tests." ON)

option(sr-index_build_tools "Build all of sr-index's tools." ON)

option(sr-index_build_benchmarks "Build all of sr-index's benchmarks." ON)


cmake_minimum_required(VERSION 3.0)


#Project Information
project(sr-index)


# Location of internal cmake modules
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${PROJECT_SOURCE_DIR}/cmake)

# Guard against in-source builds and bad build-type strings
include(ConfigSafeGuards)


#Global Setup
set(CMAKE_CXX_STANDARD 17)


# Set common include folder for module
set(COMMON_INCLUDES
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_INSTALL_PREFIX}/include
        ${CMAKE_PREFIX_PATH}/include)

include_directories(${COMMON_INCLUDES})

#TODO Verify these compilation flags
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

set(SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/include/sr-index/config.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/alphabet.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/lf.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/bwt.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/sequence_ops.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/phi.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/tools.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/definitions.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/huff_string.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/sparse_sd_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/sparse_hyb_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/rle_string.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/psi.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/io.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/enc_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/coder.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/toehold.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/sampling.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/construct_base.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/construct_sdsl.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/construct_big_bwt.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/construct.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/index_base.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/r_index.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/sr_index.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/r_csa.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/sr_csa.h
        )

find_package(sdsl)
if (NOT sdsl_FOUND)
    include(ConfigSDSL)
endif ()
include(ConfigBigBWT)
add_library(sr-index INTERFACE)
target_sources(sr-index INTERFACE ${SOURCE_FILES})
target_include_directories(sr-index INTERFACE ${PROJECT_SOURCE_DIR}/include/sr-index)
add_definitions(-DBIGBWT_EXE=\"${BigBWT_BINARY_DIR}/bigbwt\")

include(internal_utils)


if (sr-index_build_tests)
    enable_testing()

    include(ConfigGTest)

    cxx_test_with_flags_and_args(lf_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/lf_tests.cpp)
    cxx_test_with_flags_and_args(bwt_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/bwt_tests.cpp)
    cxx_test_with_flags_and_args(sequence_ops_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/sequence_ops_tests.cpp)
    cxx_test_with_flags_and_args(rle_string_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/rle_string_tests.cpp)
    cxx_test_with_flags_and_args(phi_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/phi_tests.cpp)
    cxx_test_with_flags_and_args(psi_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/psi_tests.cpp)
    cxx_test_with_flags_and_args(toehold_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/toehold_tests.cpp)
    cxx_test_with_flags_and_args(sampling_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/sampling_test.cpp)
    cxx_test_with_flags_and_args(locate_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/locate_tests.cpp)
endif ()


if (sr-index_build_tools)
endif ()


if (sr-index_build_benchmarks)
    include(ConfigGBenchmark)

    add_subdirectory(benchmark)
endif ()
