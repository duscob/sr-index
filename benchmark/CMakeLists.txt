# Configure benchmarks
set(benchmark_LIBS sr-index ${LIBS})

include(ConfigGBenchmark)
list(APPEND benchmark_LIBS benchmark::benchmark benchmark::benchmark_main)

list(APPEND benchmark_LIBS ${GFLAGS_LIB})
find_package(Threads)
list(APPEND benchmark_LIBS ${CMAKE_THREAD_LIBS_INIT})

#cxx_executable_with_flags(bm_build_ri_items "" "${benchmark_LIBS}" bm_build_ri_items.cpp)
#cxx_executable_with_flags(bm_locate "" "${benchmark_LIBS}" bm_locate.cpp)

add_subdirectory(sr-csa)
add_subdirectory(sr-index)
