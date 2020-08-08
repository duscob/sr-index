option(ri_build_tests "Build all of ri's tests." ON)

option(ri_build_tools "Build all of ri's tools." ON)

option(ri_build_benchmarks "Build all of ri's benchmarks." ON)


cmake_minimum_required(VERSION 3.0)


#Project Information
project(ri)


# Location of internal cmake modules
set(CMAKE_MODULE_PATH
        ${CMAKE_MODULE_PATH}
        ${PROJECT_SOURCE_DIR}/cmake)

# Guard against in-source builds and bad build-type strings
include(ConfigSafeGuards)


if (ri_build_tests)
    enable_testing()

    include(ConfigGTest)
endif ()


if (ri_build_tools)
endif ()


if (ri_build_benchmarks)
endif ()
