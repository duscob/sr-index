# Defines functions and macros useful for building libraries, executables
# and tests.
#
# Note:
#
# - The functions/macros defined in this file were extracted from GoogleTest
#   and adapted to our requirements.


# cxx_executable_with_flags(name cxx_flags libs srcs...)
#
# creates a named C++ executable that depends on the given libraries and
# is built from the given source files with the given compiler flags.
function(cxx_executable_with_flags name cxx_flags libs)
    add_executable(${name} ${ARGN})
    if (MSVC)
        # BigObj required for tests.
        set(cxx_flags "${cxx_flags} -bigobj")
    endif()
    if (cxx_flags)
        set_target_properties(${name}
                PROPERTIES
                COMPILE_FLAGS "${cxx_flags}")
    endif()
    if (BUILD_SHARED_LIBS)
        set_target_properties(${name}
                PROPERTIES
                COMPILE_DEFINITIONS "GTEST_LINKED_AS_SHARED_LIBRARY=1")
    endif()
    # To support mixing linking in static and dynamic libraries, link each
    # library in with an extra call to target_link_libraries.
    foreach (lib "${libs}")
        target_link_libraries(${name} ${lib})
    endforeach()
endfunction()

# cxx_executable(name dir lib srcs...)
#
# creates a named target that depends on the given libs and is built
# from the given source files.  dir/name.cc is implicitly included in
# the source file list.
function(cxx_executable name dir libs)
    cxx_executable_with_flags(
            ${name} "${cxx_default}" "${libs}" "${dir}/${name}.cc" ${ARGN})
endfunction()


# cxx_test_with_flags(name cxx_flags libs srcs...)
#
# creates a named C++ test that depends on the given libs and is built
# from the given source files with the given compiler flags.
function(cxx_test_with_flags name cxx_flags libs)
    cxx_executable_with_flags(${name} "${cxx_flags}" "${libs}" ${ARGN})
    if (WIN32 OR MINGW)
        add_test(NAME ${name}
                COMMAND "powershell" "-Command" "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/RunTest.ps1" "$<TARGET_FILE:${name}>")
    else()
        add_test(NAME ${name}
                COMMAND "$<TARGET_FILE:${name}>")
    endif()
endfunction()

# cxx_test(name libs srcs...)
#
# creates a named test target that depends on the given libs and is
# built from the given source files.  Unlike cxx_test_with_flags,
# test/name.cc is already implicitly included in the source file list.
function(cxx_test name libs)
    cxx_test_with_flags("${name}" "${cxx_default}" "${libs}"
            "test/${name}.cc" ${ARGN})
endfunction()


# cxx_test_with_flags_and_args(name cxx_flags libs cmd_args srcs...)
#
# creates a named C++ test that depends on the given libs and is built
# from the given source files with the given compiler flags.
function(cxx_test_with_flags_and_args name cxx_flags libs cmd_args)
    cxx_executable_with_flags(${name} "${cxx_flags}" "${libs}" ${ARGN})
    if (WIN32 OR MINGW)
        add_test(NAME ${name}
                COMMAND "powershell" "-Command" "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/RunTest.ps1" "$<TARGET_FILE:${name}>"
                ${cmd_args})
    else()
        add_test(NAME ${name}
                COMMAND "$<TARGET_FILE:${name}>"
                ${cmd_args})
    endif()
endfunction()
