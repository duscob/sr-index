set(ExternalProjectName json)

include(FetchContent)

FetchContent_Declare(
        ${ExternalProjectName}
        #        GIT_REPOSITORY https://github.com/nlohmann/json.git
        #        GIT_TAG v3.12.0
        URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP true
)

FetchContent_MakeAvailable(${ExternalProjectName})
#FetchContent_GetProperties(${ExternalProjectName})
#FetchContent_Populate(${ExternalProjectName})
#add_subdirectory("${${ExternalProjectName}_SOURCE_DIR}" "${${ExternalProjectName}_BINARY_DIR}")
#
#include_directories(${${ExternalProjectName}_SOURCE_DIR})
