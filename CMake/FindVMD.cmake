# Copyright (c) 2017-2019, Michael P. Howard
# This file is part of the gsd-vmd project, released under the Modified BSD License.
#
# Search for a standard VMD plugin, using ideas borrowed from FindVMD.cmake in
# GROMACS, https://github.com/gromacs/gromacs/blob/master/cmake/FindVMD.cmake
#
# The plugin directory will be searched for using the VMDDIR environment variable,
# which points to the location where the plugin libraries can be found, for
# example, /usr/local/lib/vmd on Linux. The search priority is:
#   -DVMDDIR : user-defined VMDDIR on the command line
#   $VMDDIR  : environment variable VMDDIR
#   The vmd directory pointed to by the VMD executable (this may only work on windows)
#   Standard search paths for the OS
#
# This module provides the following variables:
#   VMD_FOUND - if VMD is found
#   VMD_PLUGIN_INCLUDE_PATH - path to plugin include
#   VMD_PLUGIN_MOLFILE_PATH - path to molfile libraries

if (VMD_FOUND)
    return()
endif (VMD_FOUND)

# default guesses for where to find VMD, lowest priority
if (LINUX)
    set(VMD_PATHS "/usr/local/lib/vmd/plugins/*/molfile")
elseif (APPLE)
    set(VMD_PATHS "/Applications/*/Contents/vmd/plugins/*/molfile")
else (LINUX)
    # I don't know where VMD installs on Windows, so we don't make a default path
    set(VMD_PATHS "")
endif (LINUX)

# try to use an executable to get the vmd directory
find_program(VMD_EXECUTABLE NAMES vmd PATH_SUFFIXES bin DOC "VMD executable")
if (VMD_EXECUTABLE)
    file(STRINGS "${VMD_EXECUTABLE}" VMD_EXEC_DIR REGEX "^defaultvmddir=.*$")
    string(REGEX REPLACE "(^.*=\"?|\"$)" "" VMD_EXEC_DIR "${VMD_EXEC_DIR}")
    list(INSERT VMD_PATHS 0 "${VMD_EXEC_DIR}/plugins/*/molfile")
    unset(VMD_EXEC_DIR)
endif (VMD_EXECUTABLE)
unset(VMD_EXECUTABLE CACHE)

# use the environment variable
if (DEFINED ENV{VMDDIR})
    list(INSERT VMD_PATHS 0 "$ENV{VMDDIR}/plugins/*/molfile")
endif (DEFINED ENV{VMDDIR})

# command line definition is highest priority
if (DEFINED VMDDIR)
    list(INSERT VMD_PATHS 0 "${VMDDIR}/plugins/*/molfile")
endif (DEFINED VMDDIR)

# Attempt to find a standard VMD plugin, which will give us the plugin installation path
find_library(VMD_XYZ_PLUGIN NAME "xyzplugin.so" PATHS ${VMD_PATHS})
if (VMD_XYZ_PLUGIN)
    get_filename_component(VMD_PLUGIN_MOLFILE_PATH ${VMD_XYZ_PLUGIN} PATH)
    set(VMD_PLUGIN_MOLFILE_PATH ${VMD_PLUGIN_MOLFILE_PATH} CACHE PATH "Path to VMD molfile plugin libraries.")
endif (VMD_XYZ_PLUGIN)
unset(VMD_XYZ_PLUGIN CACHE)

# Try to get the include path if the plugin path is set
if (VMD_PLUGIN_MOLFILE_PATH)
    # the plugin path should be two directories above the molfile library path
    get_filename_component(VMD_PLUGIN_PATH ${VMD_PLUGIN_MOLFILE_PATH} PATH)
    get_filename_component(VMD_PLUGIN_PATH ${VMD_PLUGIN_PATH} PATH)

    # validate that directory exists and is real
    if (EXISTS ${VMD_PLUGIN_PATH}/include)
        set(VMD_PLUGIN_INCLUDE_PATH "${VMD_PLUGIN_PATH}/include" CACHE PATH "Path to VMD plugin includes.")
    endif (EXISTS ${VMD_PLUGIN_PATH}/include)
endif (VMD_PLUGIN_MOLFILE_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VMD DEFAULT_MSG VMD_PLUGIN_INCLUDE_PATH VMD_PLUGIN_MOLFILE_PATH)
