
find_library(GFLAGS_LIB gflags REQUIRED)
find_package(Threads)

cxx_executable_with_flags(bm_build_ri_items "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_build_ri_items.cpp)
cxx_executable_with_flags(bm_locate "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_locate.cpp)

add_subdirectory(sr-csa)
add_subdirectory(sr-index)
