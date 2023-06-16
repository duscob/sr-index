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
download_project(PROJ BigBWT
        GIT_REPOSITORY https://github.com/duscob/Big-BWT.git
        GIT_TAG master
        ${UPDATE_DISCONNECTED_IF_AVAILABLE})

add_subdirectory(${BigBWT_SOURCE_DIR} ${BigBWT_BINARY_DIR})

include_directories("${BigBWT_SOURCE_DIR}")
