
cxx_executable_with_flags(bm_construct_csa "" "${benchmark_LIBS}" bm_construct_csa.cpp)
cxx_executable_with_flags(bm_locate_csa "" "${benchmark_LIBS}" bm_locate_csa.cpp factory.h)
cxx_executable_with_flags(bm_count_csa "" "${benchmark_LIBS}" bm_count_csa.cpp factory.h)
