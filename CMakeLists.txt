option(ri_build_tests "Build all of ri's tests." ON)

option(ri_build_tools "Build all of ri's tools." ON)

option(ri_build_benchmarks "Build all of ri's benchmarks." ON)


cmake_minimum_required(VERSION 3.0)


#Project Information
project(ri)


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

set(SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/include/ri/bwt.h
        ${PROJECT_SOURCE_DIR}/include/ri/definitions.hpp
        ${PROJECT_SOURCE_DIR}/include/ri/huff_string.hpp
        ${PROJECT_SOURCE_DIR}/include/ri/sparse_sd_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/ri/sparse_hyb_vector.hpp
        ${PROJECT_SOURCE_DIR}/include/ri/rle_string.hpp)

add_library(ri INTERFACE)
target_sources(ri INTERFACE ${SOURCE_FILES})
target_include_directories(ri INTERFACE ${PROJECT_SOURCE_DIR}/include/ri)

include(internal_utils)


if (ri_build_tests)
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
endif ()


if (ri_build_tools)
endif ()


if (ri_build_benchmarks)
    include(ConfigGBenchmark)

    find_library(GFLAGS_LIB gflags REQUIRED)
    find_package(Threads)

    find_library(SDSL_LIB sdsl REQUIRED)
    find_library(DIVSUFSORT_LIB divsufsort REQUIRED)
    find_library(DIVSUFSORT64_LIB divsufsort64 REQUIRED)
    set(SDSL_LIBS ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})
    set(LIBS ${SDSL_LIBS})

    cxx_executable_with_flags(bm_build_ri_items "" "${GFLAGS_LIB};benchmark;ri;${LIBS};${CMAKE_THREAD_LIBS_INIT}" benchmark/bm_build_ri_items.cpp)
endif ()
