// Copyright (c) 2017-2020, Michael P. Howard
// Copyright (c) 2021-2026, Auburn University
// Part of gsd-vmd, released under the BSD 3-Clause License.

/*!
 * \file test_dynamic_library.cc
 * \brief Test DynamicLibrary object used for testing plugin
 */

#include "DynamicLibrary.h"
#include "catch.hpp"

//! Test the dynamic library loader for the testing harness
TEST_CASE("DynamicLibrary")
    {
    // valid plugin
    SECTION("gsdplugin")
        {
        DynamicLibrary lib("gsdplugin.so");

        // open the plugin library
        REQUIRE(lib.open());

        // initialize the plugin library
        void* ifunc = lib.load("vmdplugin_init");
        REQUIRE(ifunc);

        void* rfunc = lib.load("vmdplugin_register");
        REQUIRE(rfunc);

        void* ffunc = lib.load("vmdplugin_fini");
        REQUIRE(ffunc);

        lib.close();
        }

    // invalid plugin
    SECTION("foobar")
        {
        DynamicLibrary lib("foobar.so");

        // open the plugin library
        REQUIRE(!lib.open());

        // initialize the plugin library
        void* ifunc = lib.load("vmdplugin_init");
        REQUIRE(!ifunc);

        void* rfunc = lib.load("vmdplugin_register");
        REQUIRE(!rfunc);

        void* ffunc = lib.load("vmdplugin_fini");
        REQUIRE(!ffunc);

        lib.close();
        }
    }
