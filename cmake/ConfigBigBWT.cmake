set(ExternalProjectName BigBWT)

include(ExternalProject)
ExternalProject_Add(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/Big-BWT.git
        GIT_TAG master
        UPDATE_DISCONNECTED ON
        CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        INSTALL_COMMAND ""
)

ExternalProject_Get_property(${ExternalProjectName} BINARY_DIR)
set(${ExternalProjectName}_BINARY_DIR ${BINARY_DIR})
ExternalProject_Get_property(${ExternalProjectName} SOURCE_DIR)
include_directories(${${ExternalProjectName}_SOURCE_DIR})
