// Copyright (c) 2017-2020, Michael P. Howard
// This file is part of the gsd-vmd project, released under the Modified BSD License.

#include "GSDPlugin.h"
#include "molfile_plugin.h"

int main(int argc, char *argv[])
    {
    GSDPlugin p;
    molfile_plugin_t *plugin = p();
    if (!plugin) return 1;

    int natoms = 0;
    void *v = plugin->open_file_read("test.gsd", "gsd", &natoms);
    if (!v) return 1;

    // structure
        {
        molfile_atom_t *atoms = new molfile_atom_t[natoms];
        int optflags = 0;
        int retval = plugin->read_structure(v, &optflags, atoms);
        delete[] atoms;
        }

    // bonds
        {
        int nbonds, nbondtypes;
        int *from, *to, *bondtype;
        float *bondorder;
        char **bondtypename;
        plugin->read_bonds(v, &nbonds, &from, &to, &bondorder, &bondtype, &nbondtypes, &bondtypename);
        }

    // metadata
        {
        molfile_timestep_metadata_t *meta = new molfile_timestep_metadata_t;
        plugin->read_timestep_metadata(v, meta);
        delete meta;
        }

    // timestep
        {
        molfile_timestep_t *ts = new molfile_timestep_t;
        ts->coords = new float[natoms*3];
        ts->velocities = new float[natoms*3];

        // no more REQUIRE until memory is freed
        plugin->read_next_timestep(v, natoms, ts);

        delete[] ts->coords;
        delete[] ts->velocities;
        delete ts;
        }
    plugin->close_file_read(v);

    return 0;
    }
