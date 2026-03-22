set(ExternalProjectName sdsl)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/duscob/sdsl-lite.git
        GIT_TAG feature/io-cache
        FIND_PACKAGE_ARGS
)

set(SDSL_ENABLE_TESTS OFF CACHE BOOL "SDSL_ENABLE_TESTS")
set(SDSL_ENABLE_TUTORIALS OFF CACHE BOOL "SDSL_ENABLE_TUTORIALS")
set(SDSL_ENABLE_EXAMPLES OFF CACHE BOOL "SDSL_ENABLE_EXAMPLES")

FetchContent_MakeAvailable(${ExternalProjectName})

FetchContent_GetProperties(${ExternalProjectName})
include_directories(${${ExternalProjectName}_SOURCE_DIR}/include)
