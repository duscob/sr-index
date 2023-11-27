set(ExternalProjectName json)

include(FetchContent)

FetchContent_Declare(
        ${ExternalProjectName}
        URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
        DOWNLOAD_EXTRACT_TIMESTAMP true
)

FetchContent_MakeAvailable(${ExternalProjectName})

include_directories("${${ExternalProjectName}_SOURCE_DIR}")
