# Adapted from https://github.com/Crascit/DownloadProject/blob/master/CMakeLists.txt
#
# CAVEAT: use DownloadProject.cmake
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (CMAKE_VERSION VERSION_LESS 3.2)
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
else ()
    set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
endif ()

include(DownloadProject)
download_project(PROJ sdsl
        GIT_REPOSITORY https://github.com/duscob/sdsl-lite.git
        GIT_TAG master
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})

set(SDSL_ENABLE_TESTS OFF CACHE BOOL "SDSL_ENABLE_TESTS")
set(SDSL_ENABLE_TUTORIALS OFF CACHE BOOL "SDSL_ENABLE_TUTORIALS")
set(SDSL_ENABLE_EXAMPLES OFF CACHE BOOL "SDSL_ENABLE_EXAMPLES")

add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)

include_directories("${sdsl_SOURCE_DIR}/include")
