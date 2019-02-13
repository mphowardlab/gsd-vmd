// Copyright (c) 2017-2019, Michael P. Howard
// This file is part of the gsd-vmd project, released under the Modified BSD License.

#include "catch.hpp"
#include "GSDPlugin.h"
#include "molfile_plugin.h"

#include <string>

//! Test the GSD plugin
TEST_CASE("GSDPlugin")
    {
    GSDPlugin p;
    molfile_plugin_t *plugin = p();
    REQUIRE(plugin);

    CHECK(plugin->abiversion == vmdplugin_ABIVERSION);
    CHECK(plugin->type == std::string(MOLFILE_PLUGIN_TYPE));
    CHECK(std::string(plugin->name) == "gsd");
    CHECK(std::string(plugin->prettyname) == "HOOMD-blue GSD File");
    CHECK(std::string(plugin->author) == "Michael P. Howard");
    CHECK(plugin->majorv == 0);
    CHECK(plugin->minorv == 1);
    CHECK(plugin->is_reentrant == VMDPLUGIN_THREADSAFE);
    CHECK(std::string(plugin->filename_extension) == "gsd");

    CHECK(plugin->open_file_read);
    CHECK(plugin->read_structure);
    CHECK(plugin->read_bonds);
    CHECK(plugin->read_next_timestep);
    CHECK(plugin->read_timestep_metadata);
    CHECK(plugin->close_file_read);
    }

//! Test reading a representative GSD file
TEST_CASE("Read GSD")
    {
    GSDPlugin p;
    molfile_plugin_t *plugin = p();
    REQUIRE(plugin);

    int natoms = 0;
    void *v = plugin->open_file_read("test.gsd", "gsd", &natoms);
    REQUIRE(v != NULL);
    REQUIRE(natoms == 2);

    SECTION("structure")
        {
        molfile_atom_t *atoms = new molfile_atom_t[natoms];
        REQUIRE(atoms);

        // cannot use any REQUIRE statements until atoms has been freed
        int optflags = 0;
        int retval = plugin->read_structure(v, &optflags, atoms);
        CHECK(retval == MOLFILE_SUCCESS);
        CHECK( (optflags & MOLFILE_MASS) );
        CHECK( (optflags & MOLFILE_CHARGE) );
        CHECK( (optflags & MOLFILE_RADIUS) );

        // mass, charge, radius of particle 2
        CHECK(std::string(atoms[0].type) == "A");
        CHECK(std::string(atoms[0].name) == "A");
        CHECK(atoms[0].mass == Approx(6.));
        CHECK(atoms[0].charge == Approx(-1.));
        CHECK(atoms[0].radius == Approx(1.));
        // resname, resid, segid, and chain should all be unset
        CHECK(std::string(atoms[0].resname) == std::string());
        CHECK(atoms[0].resid == 0);
        CHECK(std::string(atoms[0].segid) == std::string());
        CHECK(std::string(atoms[0].chain) == std::string());

        // mass, charge, radius of particle 2
        CHECK(std::string(atoms[1].type) == "B");
        CHECK(std::string(atoms[1].name) == "B");
        CHECK(atoms[1].mass == Approx(8.));
        CHECK(atoms[1].charge == Approx(1.));
        CHECK(atoms[1].radius == Approx(2.));
        // resname, resid, segid, and chain should all be unset
        CHECK(std::string(atoms[1].resname) == std::string());
        CHECK(atoms[1].resid == 0);
        CHECK(std::string(atoms[1].segid) == std::string());
        CHECK(std::string(atoms[1].chain) == std::string());

        delete[] atoms;
        }

    SECTION("bonds")
        {
        int nbonds, nbondtypes;
        int *from, *to, *bondtype;
        float *bondorder;
        char **bondtypename;
        int retval = plugin->read_bonds(v, &nbonds, &from, &to, &bondorder, &bondtype, &nbondtypes, &bondtypename);
        REQUIRE(retval == MOLFILE_SUCCESS);

        // bondorder is not supplied
        REQUIRE(bondorder == NULL);

        // bond from 0->1, which is 1-indexed in vmd
        REQUIRE(nbonds == 1);
        CHECK(from[0] == 1);
        CHECK(to[0] == 2);

        // only one bond type
        REQUIRE(nbondtypes == 1);
        CHECK(std::string(bondtypename[0]) == "C");
        }

    SECTION("metadata")
        {
        molfile_timestep_metadata_t *meta = new molfile_timestep_metadata_t;

        // no more REQUIRE until meta is freed
        int retval = plugin->read_timestep_metadata(v, meta);
        CHECK(retval == MOLFILE_SUCCESS);
        CHECK(meta->count == 1);
        CHECK(meta->has_velocities);

        delete meta;
        }

    SECTION("timestep")
        {
        molfile_timestep_t *ts = new molfile_timestep_t;
        ts->coords = new float[natoms*3];
        ts->velocities = new float[natoms*3];

        // no more REQUIRE until memory is freed
        int retval = plugin->read_next_timestep(v, natoms, ts);
        CHECK(retval == MOLFILE_SUCCESS);

        // check timestep
        CHECK(ts->physical_time == Approx(0.));

        // check simulation box
        CHECK(ts->A == Approx(20.0));
        CHECK(ts->B == Approx(20.0));
        CHECK(ts->C == Approx(20.0));
        CHECK(ts->alpha == Approx(90.0));
        CHECK(ts->beta == Approx(90.0));
        CHECK(ts->gamma == Approx(90.0));

        // check positions
        // particle 1: x, y, z
        CHECK(ts->coords[0] == Approx(1.));
        CHECK(ts->coords[1] == Approx(2.));
        CHECK(ts->coords[2] == Approx(3.));
        // particle 2: x, y, z
        CHECK(ts->coords[3] == Approx(4.));
        CHECK(ts->coords[4] == Approx(5.));
        CHECK(ts->coords[5] == Approx(6.));

        // check velocities
        // particle 1: x, y, z
        CHECK(ts->velocities[0] == Approx(-1.));
        CHECK(ts->velocities[1] == Approx(-2.));
        CHECK(ts->velocities[2] == Approx(-3.));
        // particle 2: x, y, z
        CHECK(ts->velocities[3] == Approx(-4.));
        CHECK(ts->velocities[4] == Approx(-5.));
        CHECK(ts->velocities[5] == Approx(-6.));

        // passing a null pointer should return success, but we've also reached EOF
        // so it should just return as such
        retval = plugin->read_next_timestep(v, natoms, NULL);
        CHECK(retval == MOLFILE_EOF);

        delete[] ts->coords;
        delete[] ts->velocities;
        delete ts;
        }

    plugin->close_file_read(v);
    }
