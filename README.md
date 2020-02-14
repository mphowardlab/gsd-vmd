# HOOMD-blue GSD plugin for VMD

Provides a VMD molfile plugin reader for GSD files generated by HOOMD-blue.
In order to compile the plugin, you will first need a working installation
of VMD. Then, installation can be as simple as

```bash
cd /path/to/gsd-vmd
mkdir build && cd build
cmake ..
make install
```

On Linux, gsd-vmd should automatically detect the location of your VMD
installation, provided that it is on your path. On macOS, the VMD
installation can be trickier to find. You can hint the location
of your VMD plugin directory by setting `VMDDIR` either as an
environment variable or as a build-time definition. Failing this,
you can also force the values of your plugin installation directory
and plugin include directory with the `VMD_PLUGIN_INCLUDE_PATH`
and `VMD_PLUGIN_MOLFILE_PATH`. These directories are the location
of your plugin headers (e.g., `molfile_plugin.h`) and
system-specific molfile libraries, respectively.

By default, gsd-vmd will be installed into `VMD_PLUGIN_MOLFILE_PATH`.
You can specify an alternative installation location; however, you
must ensure that this is added to your VMD search path, which is a
tricky endeavour and not recommended. To uninstall the library, simply
remove `gsdplugin.so` from the installation location.

This plugin has been tested on Linux using gcc and clang compilers
for x86 architecture. It has also been successfully built for versions
of the Mac operating system supporting 32-bit applications. However,
beginning with macOS 10.15, VMD is currently not officially supported
on macOS because 32-bit applications were deprecated. Users have
reported success building the plugin against an unofficial 64-bit port
of VMD. If you would like to try this, you should set the CMake
option `CMAKE_OSX_ARCHITECTURES=x86_64` to only build in 64-bit.
No support has been tested for Windows.

To test your build (on Linux), run

```bash
make test
```

out of your build directory. If you have valgrind installed, a memcheck
will also be run to detect memory leaks.

## Requirements
* A compiler supporting basic C99 standard (tested gcc and clang)
* CMake >= 3.1
* valgrind (optional, for testing only)

## Source code
**GSD**: The GSD library (https://github.com/glotzerlab/gsd) is used for
file reading under the following license:

    Copyright (c) 2016-2020 The Regents of the University of Michigan
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
    ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**Catch**: The Catch framework is used for unit test validation under
the Boost Software License (version 1.0).
