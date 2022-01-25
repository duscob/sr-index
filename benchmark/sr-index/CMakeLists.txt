
cxx_executable_with_flags(bm_construct_ri "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_construct_ri.cpp)
cxx_executable_with_flags(bm_locate_ri "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_locate_ri.cpp factory.h)
