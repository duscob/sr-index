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
set(CMAKE_CXX_STANDARD 14)


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
        ${PROJECT_SOURCE_DIR}/include/sr-index/bwt.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/predecessor.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/phi.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/tools.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/r_index.h
        ${PROJECT_SOURCE_DIR}/include/sr-index/definitions.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/huff_string.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/sparse_sd_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/sparse_hyb_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/sr-index/rle_string.hpp)

add_library(sr-index INTERFACE)
target_sources(sr-index INTERFACE ${SOURCE_FILES})
target_include_directories(sr-index INTERFACE ${PROJECT_SOURCE_DIR}/include/sr-index)

include(internal_utils)


if (sr-index_build_tests)
    enable_testing()

    include(ConfigGTest)

    find_library(SDSL_LIB sdsl REQUIRED)
    find_library(DIVSUFSORT_LIB divsufsort REQUIRED)
    find_library(DIVSUFSORT64_LIB divsufsort64 REQUIRED)
    set(SDSL_LIBS ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})
    set(LIBS ${SDSL_LIBS})

    cxx_test_with_flags_and_args(bwt_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/bwt_tests.cpp)
    cxx_test_with_flags_and_args(predecessor_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/predecessor_tests.cpp)
    cxx_test_with_flags_and_args(rle_string_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/rle_string_tests.cpp)
    cxx_test_with_flags_and_args(phi_tests "" "gtest;gtest_main;gmock;${LIBS}" "" ${PROJECT_SOURCE_DIR}/test/phi_tests.cpp)
endif ()


if (sr-index_build_tools)
endif ()


if (sr-index_build_benchmarks)
    include(ConfigGBenchmark)

    find_library(GFLAGS_LIB gflags REQUIRED)
    find_package(Threads)

    find_library(SDSL_LIB sdsl REQUIRED)
    find_library(DIVSUFSORT_LIB divsufsort REQUIRED)
    find_library(DIVSUFSORT64_LIB divsufsort64 REQUIRED)
    set(SDSL_LIBS ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})
    set(LIBS ${SDSL_LIBS})

    cxx_executable_with_flags(bm_build_ri_items "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" benchmark/bm_build_ri_items.cpp)
    cxx_executable_with_flags(bm_locate "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" benchmark/bm_locate.cpp)
endif ()
