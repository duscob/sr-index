option(ri_build_tests "Build all of ri's tests." ON)

option(ri_build_tools "Build all of ri's tools." ON)

option(ri_build_benchmarks "Build all of ri's benchmarks." ON)


cmake_minimum_required(VERSION 3.0)


#Project Information
project(ri)


if (ri_build_tests)
    enable_testing()
endif ()


if (ri_build_tools)
endif ()


if (ri_build_benchmarks)
endif ()
