
cxx_executable_with_flags(bm_construct_csa "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_construct_csa.cpp index_base.h csa.h sr_csa.h csa_raw.h)
cxx_executable_with_flags(bm_locate_csa "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_locate_csa.cpp csa.h factory.h)
