# Configure benchmarks
set(benchmark_LIBS sr-index ${LIBS})

include(ConfigGBenchmark)
list(APPEND benchmark_LIBS benchmark::benchmark benchmark::benchmark_main)

list(APPEND benchmark_LIBS ${GFLAGS_LIB})
find_package(Threads)
list(APPEND benchmark_LIBS ${CMAKE_THREAD_LIBS_INIT})


add_subdirectory(sr-csa)
add_subdirectory(sr-index)
