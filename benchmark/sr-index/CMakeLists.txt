
cxx_executable_with_flags(bm_construct_ri "" "${benchmark_LIBS}" bm_construct_ri.cpp)
cxx_executable_with_flags(bm_locate_ri "" "${benchmark_LIBS}" bm_locate_ri.cpp factory.h)
cxx_executable_with_flags(bm_count_ri "" "${benchmark_LIBS}" bm_count_ri.cpp factory.h)
