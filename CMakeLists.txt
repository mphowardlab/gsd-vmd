# Copyright (c) 2017-2020, Michael P. Howard
# Copyright (c) 2021, Auburn University
# This file is part of the gsd-vmd project, released under the Modified BSD License.

cmake_minimum_required(VERSION 3.10)

project(gsd-vmd)
set(CMAKE_C_STANDARD 99)

# find the vmd directory, and configure installation prefix
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)
find_package(VMD REQUIRED)
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${VMD_PLUGIN_MOLFILE_PATH} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# configure library target
include_directories(${VMD_PLUGIN_INCLUDE_PATH})
add_library(gsdplugin SHARED gsdplugin.c gsd.c)
set_target_properties(gsdplugin PROPERTIES PREFIX "" SUFFIX ".so")
install(TARGETS gsdplugin DESTINATION ${VMD_PLUGIN_MOLFILE_PATH})

# add the testing tree
option(BUILD_TESTING "Build the testing tree." OFF)
enable_testing()
include(CTest)
if (BUILD_TESTING)
    add_subdirectory(tests)
endif (BUILD_TESTING)
