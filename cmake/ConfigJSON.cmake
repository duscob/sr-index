set(ExternalProjectName json)

include(FetchContent)

FetchContent_Declare(
        ${ExternalProjectName}
        DOWNLOAD_EXTRACT_TIMESTAMP true
        URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)

#FetchContent_MakeAvailable(${ExternalProjectName})
FetchContent_GetProperties(${ExternalProjectName})
FetchContent_Populate(${ExternalProjectName})
add_subdirectory("${${ExternalProjectName}_SOURCE_DIR}" "${${ExternalProjectName}_BINARY_DIR}")

include_directories(${${ExternalProjectName}_SOURCE_DIR})
