
find_library(GFLAGS_LIB gflags REQUIRED)
find_package(Threads)

find_library(SDSL_LIB sdsl REQUIRED)
find_library(DIVSUFSORT_LIB divsufsort REQUIRED)
find_library(DIVSUFSORT64_LIB divsufsort64 REQUIRED)
set(SDSL_LIBS ${SDSL_LIB} ${DIVSUFSORT_LIB} ${DIVSUFSORT64_LIB})
set(LIBS ${SDSL_LIBS})

cxx_executable_with_flags(bm_build_ri_items "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_build_ri_items.cpp)
cxx_executable_with_flags(bm_locate "" "${GFLAGS_LIB};benchmark;sr-index;${LIBS};${CMAKE_THREAD_LIBS_INIT}" bm_locate.cpp)

add_subdirectory(sr-csa)