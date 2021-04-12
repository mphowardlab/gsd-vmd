# Copyright (c) 2017-2020, Michael P. Howard
# Copyright (c) 2021, Auburn University
# This file is part of the gsd-vmd project, released under the Modified BSD License.

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

add_library(catch catch.cc)

add_library(gsdplugin_loader DynamicLibrary.cc GSDPlugin.cc)
target_link_libraries(gsdplugin_loader dl)

add_executable(test_gsdplugin test_gsdplugin.cc test_dynamic_library.cc)
target_link_libraries(test_gsdplugin gsdplugin gsdplugin_loader catch)
add_test(gsdplugin test_gsdplugin)

find_program(VALGRIND valgrind)
if(VALGRIND)
    add_executable(memcheck_gsdplugin memcheck_gsdplugin.cc)
    target_link_libraries(memcheck_gsdplugin gsdplugin gsdplugin_loader)
    add_test(NAME memcheck-gsdplugin COMMAND ${VALGRIND} --leak-check=full --error-exitcode=123 ${CMAKE_CURRENT_BINARY_DIR}/memcheck_gsdplugin)
endif()

add_custom_command(TARGET test_gsdplugin POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                           ${CMAKE_CURRENT_SOURCE_DIR}/test.gsd
                           ${CMAKE_CURRENT_BINARY_DIR}/test.gsd)
