set(ExternalProjectName googlebenchmark)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.0
        FIND_PACKAGE_ARGS
)

#Benchmark
# If you want to self-test benchmark lib too, turn me ON
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "BENCHMARK_ENABLE_TESTING")
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "BENCHMARK_ENABLE_GTEST_TESTS")

# Prevent Benchmark installation
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(${ExternalProjectName})
