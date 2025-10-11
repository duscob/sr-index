set(ExternalProjectName googletest)

include(FetchContent)
FetchContent_Declare(
        ${ExternalProjectName}
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.11.0
        FIND_PACKAGE_ARGS NAMES GTest
)

FetchContent_MakeAvailable(${ExternalProjectName})
