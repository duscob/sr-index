set(ExternalProjectName gflags)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/gflags/gflags.git
        GIT_TAG v2.3.0
        #        GIT_TAG main
        FIND_PACKAGE_ARGS
)

set(GFLAGS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(GFLAGS_BUILD_PACKAGING OFF CACHE BOOL "" FORCE)
set(GFLAGS_BUILD_PACKAGING OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(${ExternalProjectName})
