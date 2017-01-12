// Copyright (c) 2017, Michael P. Howard
// This file is part of the gsd-vmd project, released under the Modified BSD License.

/*!
 * \file gsdplugin.c
 * \brief VMD molfile plugin to read HOOMD-blue GSD files
 */

#include "gsd.h"
#include "molfile_plugin.h"
#include "vmdconio.h"

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

//! Macro to safely free and NULL a pointer \a p
#define SAFE_FREE(p) do{ if(p){ free(p); p = NULL; } } while(0)

//! GSD handle object
typedef struct gsd_handle gsd_handle_t;

//! String type map
/*!
 * Maps an index to a string name
 */
typedef struct
    {
    int ntypes;     //!< Number of types mapped
    char **type;    //!< Names of types
    } typemap_t;

//! Reallocate a type map and null entry strings
/*!
 * \param typemap String type map to reallocate
 * \param ntypes Number of types to allocate for
 *
 * When \a ntypes is 0, all memory is freed and pointers are nulled. When \a ntypes
 * is greater than 0, memory for the map is allocated and each entry initialized
 * to a NULL pointer.
 *
 * \warning Any values stored in the map are lost on reallocation.
 */
static void reallocate_typemap(typemap_t* typemap, int ntypes)
    {
    if (!typemap) return;

    // free any existing memory
    if(typemap->type)
        {
        for (int i=0; i < typemap->ntypes; ++i)
            {
            SAFE_FREE(typemap->type[i]);
            }
        SAFE_FREE(typemap->type);
        }

    typemap->ntypes = ntypes;
    if (ntypes > 0)
        {
        typemap->type = (char**)malloc(ntypes * sizeof(char*));
        for (int i=0; i < typemap->ntypes; ++i)
            {
            (typemap->type)[i] = NULL;
            }
        }
    }

//! Type map constructor
/*!
 * \param ntypes Number of types to hold in the map
 * \returns An initialized type map with all entries set to NULL pointers
 */
static typemap_t* allocate_typemap(int ntypes)
    {
    typemap_t *typemap = (typemap_t *)malloc(sizeof(typemap_t));

    typemap->ntypes = 0;
    typemap->type = NULL;
    reallocate_typemap(typemap, ntypes);

    return typemap;
    }

//! Type map destructor
/*!
 * \param typemap String type map to free
 *
 * \post All memory is freed, and \a typemap is set to NULL.
 *
 * This function is safe to call even if \a typemap is not allocated.
 */
static void free_typemap(typemap_t* typemap)
    {
    if (typemap)
        {
        // free any existing memory by calling a reallocation
        reallocate_typemap(typemap, 0);
        SAFE_FREE(typemap);
        }
    }

//! GSD trajectory
typedef struct
    {
    gsd_handle_t* handle;       //!< GSD file handle
    int frame;                  //!< Current frame index
    int numframes;              //!< Number of frames in gsd file
    int natoms;                 //!< Number of atoms in first frame
    typemap_t *typemap;         //!< Type map

    int nbonds;                 //!< Number of bonds
    int *bond_from;             //!< First particle in bond (1-indexed)
    int *bond_to;               //!< Second particle in bond (2-indexed)
    typemap_t *bondmap;         //!< Bond map
    } gsd_trajectory_t;

//! Constructor for GSD trajectory
/*!
 * \returns An allocated and default initialized GSD trajectory
 */
static gsd_trajectory_t* allocate_gsd_trajectory()
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t *)calloc(1,sizeof(gsd_trajectory_t));

    if (gsd)
        {
        gsd->handle = (gsd_handle_t*)malloc(sizeof(gsd_handle_t));

        gsd->frame = 0;
        gsd->numframes = 0;

        gsd->natoms = 0;
        gsd->typemap = allocate_typemap(0);

        gsd->nbonds = 0;
        gsd->bond_from = NULL;
        gsd->bond_to = NULL;
        gsd->bondmap = allocate_typemap(0);
        }

    return gsd;
    }

//! Destructor for GSD trajectory
/*!
 * \param gsd GSD trajectory object
 *
 * \post All memory is freed and \a gsd is set to NULL.
 *
 * The GSD file handle is closed, internal pointers to memory are freed, and
 * the trajectory memory is deallocated.
 *
 * This function is safe to call even if \a gsd is not allocated.
 */
static void free_gsd_trajectory(gsd_trajectory_t *gsd)
    {
    if (gsd)
        {
        if (gsd->handle)
            {
            gsd_close(gsd->handle);
            SAFE_FREE(gsd->handle);
            }
        free_typemap(gsd->typemap);

        gsd->nbonds = 0;
        if (gsd->bond_from)
            {
            SAFE_FREE(gsd->bond_from);
            }
        if (gsd->bond_to)
            {
            SAFE_FREE(gsd->bond_to);
            }
        free_typemap(gsd->bondmap);

        SAFE_FREE(gsd);
        }
    }

//! Read a chunk from the GSD file
/*!
 * \param handle Pointer to the GSD file handle
 * \param data Array to write the data into
 * \param frame Frame index
 * \param name Name of chunk to read
 * \param expected_size Bytes allocated to \a data hold the chunk
 * \param expected_N Expected number of particles
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 *
 * If \a expected_N is nonzero, then the chunk size is validated to ensure it
 * contains \a expected_N entries.
 */
static int read_chunk(gsd_handle_t *handle,
                      void *data,
                      uint64_t frame,
                      const char *name,
                      size_t expected_size,
                      unsigned int expected_N)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, frame, name);
    if (entry == NULL)
        {
        // silently ignore missing entries
        return -3;
        }
    else if (expected_N != 0 && entry->N != expected_N)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Incorrect number of entries in chunk '%s' at frame %d.\n"
                                    "           Expected %d and found %d\n", name, frame, expected_N, entry->N);
        return -2;
        }

    // check that size of entry matches the expected size
    size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
    if (actual_size != expected_size)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Incorrect entry size in chunk '%s' at frame %d.\n"
                                    "           Expected %d B and found %d B\n", name, frame, expected_size, actual_size);
        return -2;
        }

    int retval = gsd_read_chunk(handle, data, entry);
    if (retval == -1)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Chunk '%s' : %s\n", name, strerror(errno));
        }
    else if (retval != 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading chunk '%s'\n", name);
        }

    return retval;
    }

//! Open the GSD file for reading
/*!
 * \param filename GSD filename
 * \param filetype VMD supplied filetype (unused)
 * \param natoms Number of atoms in trajectory
 *
 * \returns A pointer to a GSD trajectory object with an open file handle.
 *
 * A GSD trajectory is allocated, and a file handle is safely opened. The number
 * of frames and number of particles in the GSD file are read.
 * Although GSD supports changing number of particles, VMD does not, so \a natoms
 * will be set from the value of N stored in frame 0.
 */
static void* open_gsd_read(const char *filename, const char *filetype, int *natoms)
    {
    if (!filename) return NULL;

    gsd_trajectory_t *gsd = allocate_gsd_trajectory();
    if (!gsd) return NULL;

    int retval = gsd_open(gsd->handle, filename, GSD_OPEN_READONLY);
    if (retval == -1)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) '%s': %s\n", filename, strerror(errno));
        }
    else if (retval == -2)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) '%s' is not a valid GSD file\n", filename);
        }
    else if (retval == -3)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Invalid GSD file version in '%s'\n", filename);
        }
    else if (retval == -4)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Corrupt GSD file '%s'\n", filename);
        }
    else if (retval == -5)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Unable to allocate memory opening '%s'\n", filename);
        }
    else if (retval != 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error opening '%s'\n", filename);
        }
    // safe free and return in all cases of error
    if (retval != 0)
        {
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // validate schema
    if (strcmp(gsd->handle->header.schema, "hoomd") != 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Invalid schema in '%s', expecting 'hoomd'\n", filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }
    if (gsd->handle->header.schema_version >= gsd_make_version(2,0))
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Invalid schema version in '%s', expecting 1.x\n", filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // validate that at least one frame is written (is this just given?)
    gsd->numframes = gsd_get_nframes(gsd->handle);
    if (gsd->numframes == 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) GSD file '%s' does not contain any frames\n", filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // read the number of particles
    *natoms = 0;
    read_chunk(gsd->handle, natoms, 0, "particles/N", sizeof(int), 0);
    if (*natoms == 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) No particles found in first frame of '%s'\n", filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }
    gsd->natoms = *natoms;

    return gsd;
    }

//! Read the type map from the GSD file into the trajectory
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 *
 * Reads in the type map from the GSD file as NULL terminated strings, respecting
 * the maximum character length that can be accommodated by VMD. If no type map
 * is set, the default of 0 -> A is assumed, per the HOOMD-blue schema.
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_typemap(gsd_trajectory_t *gsd, molfile_atom_t *atoms)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(gsd->handle, 0, "particles/types");
    if (entry != NULL) // types are present
        {
        size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
        char* data = (char*)malloc(actual_size);
        int retval = gsd_read_chunk(gsd->handle, data, entry);
        if (retval == -1)
            {
            SAFE_FREE(data);
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Type mapping 'particles/types' : %s\n", strerror(errno));
            return MOLFILE_ERROR;
            }
        else if (retval != 0)
            {
            SAFE_FREE(data);
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading type mapping 'particles/types'\n");
            return MOLFILE_ERROR;
            }

        // determine the maximum copy size from the molfile_atom_t
        if (atoms == NULL) return MOLFILE_ERROR;
        const size_t max_name = sizeof(atoms->name);
        const size_t max_type = sizeof(atoms->name);
        const size_t max_nametype = (max_name < max_type) ? max_name-1 : max_type-1;
        if (max_nametype < (entry->M - 1))
            {
            vmdcon_printf(VMDCON_WARN, "gsdplugin) Type names cannot exceed %d characters, truncating\n", max_nametype);
            }

        // remalloc the type mapping and copy from the gsd data with null termination
        reallocate_typemap(gsd->typemap, entry->N);
        for (int i=0; i < entry->N; ++i)
            {
            const char *name = data + i*entry->M;
            // get size of the name
            size_t l = strnlen(name, entry->M);
            if (l > max_nametype)
                {
                l = max_nametype;
                }

            (gsd->typemap->type)[i] = (char*)malloc((l+1) * sizeof(char));
            strncpy((gsd->typemap->type)[i], name, l);
            (gsd->typemap->type)[i][l] = '\0'; // force null termination
            }
        SAFE_FREE(data);
        }
    else
        {
        // initialize default types from the HOOMD spec
        reallocate_typemap(gsd->typemap, 1);
        gsd->typemap->type[0] = (char*)malloc(2*sizeof(char));
        strncpy(gsd->typemap->type[0],"A\0",2);
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle types from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_types(gsd_trajectory_t *gsd, molfile_atom_t *atoms)
    {
    if (read_gsd_typemap(gsd, atoms) != MOLFILE_SUCCESS) return MOLFILE_ERROR;
    uint32_t *typeid = (uint32_t*)calloc(gsd->natoms, sizeof(uint32_t));
    int retval = read_chunk(gsd->handle, typeid, 0, "particles/typeid", gsd->natoms * sizeof(uint32_t), gsd->natoms);
    if (retval == 0 || retval == -3)
        {
        molfile_atom_t *a = atoms;
        for (int i=0; i < gsd->natoms; ++i, ++a)
            {
            unsigned int typeid_i = typeid[i];
            strncpy(a->name, gsd->typemap->type[typeid_i], sizeof(atoms->name));
            strncpy(a->type, gsd->typemap->type[typeid_i], sizeof(atoms->type));
            }
        }
    else
        {
        SAFE_FREE(typeid);
        return MOLFILE_ERROR;
        }

    SAFE_FREE(typeid);
    return MOLFILE_SUCCESS;
    }

//! Read particle types from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 * \param tmp Pointer to temporary memory allocated to hold the number of atoms
 *            to read from \a gsd
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_mass(gsd_trajectory_t *gsd, molfile_atom_t *atoms, float *tmp)
    {
    if (!tmp) return MOLFILE_ERROR;

    int retval = read_chunk(gsd->handle, tmp, 0, "particles/mass", gsd->natoms * sizeof(float), gsd->natoms);
    if (retval == 0)
        {
        molfile_atom_t *a = atoms;
        for (int i=0; i < gsd->natoms; ++i, ++a)
            {
            a->mass = tmp[i];
            }
        }
    else if (retval != -3)
        {
        return MOLFILE_ERROR;
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle charge from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 * \param tmp Pointer to temporary memory allocated to hold the number of atoms
 *            to read from \a gsd
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_charge(gsd_trajectory_t *gsd, molfile_atom_t *atoms, float *tmp)
    {
    if (!tmp) return MOLFILE_ERROR;

    int retval = read_chunk(gsd->handle, tmp, 0, "particles/charge", gsd->natoms * sizeof(float), gsd->natoms);
    if (retval == 0)
        {
        molfile_atom_t *a = atoms;
        for (int i=0; i < gsd->natoms; ++i, ++a)
            {
            a->charge = tmp[i];
            }
        }
    else if (retval != -3)
        {
        return MOLFILE_ERROR;
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle radius from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 * \param tmp Pointer to temporary memory allocated to hold the number of atoms
 *            to read from \a gsd
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_radius(gsd_trajectory_t *gsd, molfile_atom_t *atoms, float *tmp)
    {
    if (!tmp) return MOLFILE_ERROR;

    int retval = read_chunk(gsd->handle, tmp, 0, "particles/diameter", gsd->natoms * sizeof(float), gsd->natoms);
    if (retval == 0)
        {
        molfile_atom_t *a = atoms;
        for (int i=0; i < gsd->natoms; ++i, ++a)
            {
            a->radius = 0.5f * tmp[i];
            }
        }
    else if (retval != -3)
        {
        return MOLFILE_ERROR;
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle types from the GSD file
/*!
 * \param mydata GSD trajectory
 * \param optflags VMD optional flags (output)
 * \param atoms VMD atom properties (output)
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 *
 * Atom properties are first set to their default values, per the HOOMD schema,
 * and then we attempt to read the data from frame 0.
 *
 * \sa read_gsd_types for how types are read
 * \sa read_gsd_mass for how masses are read
 * \sa read_gsd_charge for how charge is read
 * \sa read_gsd_radius for how radius is read
 */
static int read_gsd_structure(void *mydata, int *optflags, molfile_atom_t *atoms)
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t*)mydata;

    // loop through the atoms and fill in with defaults per the HOOMD GSD schema
        {
        molfile_atom_t *a = atoms;
        for (int i=0; i < gsd->natoms; ++i, ++a)
            {
            strncpy(a->name, "A\0", sizeof(atoms->name));
            strncpy(a->type, "A\0", sizeof(atoms->type));
            a->resname[0]='\0';
            a->resid=0;
            a->segid[0]='\0';
            a->chain[0]='\0';

            a->mass =    1.0f;
            a->charge =  0.0f;
            a->radius =  0.5f;
            }
        }
    *optflags = MOLFILE_MASS | MOLFILE_CHARGE | MOLFILE_RADIUS;

    // map the particle types
    if (read_gsd_types(gsd, atoms) != MOLFILE_SUCCESS) return MOLFILE_ERROR;

    // try to set optional properties
    float *props = (float*)calloc(gsd->natoms, sizeof(float));
    // mass
    if (read_gsd_mass(gsd, atoms, props) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(props);
        return MOLFILE_ERROR;
        }
    // charge
    if (read_gsd_charge(gsd, atoms, props) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(props);
        return MOLFILE_ERROR;
        }
    // radius
    if (read_gsd_radius(gsd, atoms, props) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(props);
        return MOLFILE_ERROR;
        }
    SAFE_FREE(props);

    return MOLFILE_SUCCESS;
    }

//! Read the type map from the GSD file into the trajectory
/*!
 * \param handle GSD file handle
 * \param name Name of the bond map chunk
 * \param bondmap Bond name map (output)
 *
 * \returns MOLFILE_SUCCESS or MOLFILE_ERROR on success / failure.
 *
 * Reads in the bond type map from frame 0 of the GSD file as NULL terminated strings,
 * and saves it into \a bondmap. This method may be used to read bond, angle, and dihedral names.
 */
static int read_bondmap(gsd_handle_t *handle,
                        const char *name,
                        typemap_t *bondmap)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, 0, name);
    if (entry != NULL) // types are present
        {
        size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
        char* data = (char*)malloc(actual_size);
        int retval = gsd_read_chunk(handle, data, entry);
        if (retval == -1)
            {
            SAFE_FREE(data);
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Type mapping '%s' : %s\n", name, strerror(errno));
            return MOLFILE_ERROR;
            }
        else if (retval != 0)
            {
            SAFE_FREE(data);
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading type mapping '%s'\n", name);
            return MOLFILE_ERROR;
            }

        reallocate_typemap(bondmap, entry->N);
        if (bondmap && bondmap->type)
            {
            for (int i=0; i < entry->N; ++i)
                {
                const char *name = data + i*entry->M;
                // get size of the name
                size_t l = strnlen(name, entry->M);

                // resizing guarantees that all member chars are nulled, so can just malloc
                bondmap->type[i] = (char*)malloc((l+1) * sizeof(char));
                strncpy(bondmap->type[i], name, l);
                bondmap->type[i][l] = '\0'; // force null termination
                }
            }
        SAFE_FREE(data);
        }

    return MOLFILE_SUCCESS;
    }

//! Read the bonds from the GSD file into the trajectory
/*!
 * \param mydata GSD trajectory
 * \param nbonds Pointer to VMD number of bonds (output)
 * \param from Pointer to list of first atom in bonds, 1-indexed (output)
 * \param to Point to list of second atom in bonds, 1-indexed (output)
 * \param bondorder Pointer to bond order values (output)
 * \param bondtype Pointer to list of type of each bond (output)
 * \param nbondtypes Pointer to VMD number of bond types (output)
 * \param bondtypename Pointer to VMD list of bond type names (output)
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 *
 * HOOMD GSD does not supply any bond order data, so \a bondorder is set to
 * a NULL pointer. Bonds are read in, and converted to VMD 1-indexing. Bond names
 * are read in using a typemap_t, and are stored in the GSD trajectory \a mydata.
 */
static int read_gsd_bonds(void *mydata,
                          int *nbonds,
                          int **from,
                          int **to,
                          float **bondorder,
                          int **bondtype,
                          int *nbondtypes,
                          char ***bondtypename)
    {
    // gsd does not supply a bond order
    *bondorder = NULL;

    // default is to not supply any bonds
    *nbonds = 0;
    *from = NULL;
    *to = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    // check number of bonds and exit early if no bonds are present, or on read error
    gsd_trajectory_t *gsd = (gsd_trajectory_t*)mydata;
    int retval = read_chunk(gsd->handle, &gsd->nbonds,  0, "bonds/N", sizeof(int), 0);
    if (retval == -3 || gsd->nbonds == 0)
        {
        // return successfully if bonds are not present in the file
        return MOLFILE_SUCCESS;
        }
    else if (retval != 0)
        {
        // exit with read error
        return MOLFILE_ERROR;
        }

    // acquire the bondname map
    if (read_bondmap(gsd->handle, "bonds/types", gsd->bondmap) != MOLFILE_SUCCESS)
        {
        return MOLFILE_ERROR;
        }

    // read in the bonds, and remap them with 1-indexing
    uint32_t *bonds = (uint32_t*)malloc(2*gsd->nbonds * sizeof(uint32_t));
    retval = read_chunk(gsd->handle, bonds, 0, "bonds/group", 2*gsd->nbonds*sizeof(uint32_t), gsd->nbonds);
    if (retval != 0)
        {
        SAFE_FREE(bonds);
        return MOLFILE_ERROR;
        }
    gsd->bond_from = (int*)malloc(gsd->nbonds * sizeof(int));
    gsd->bond_to = (int*)malloc(gsd->nbonds * sizeof(int));
    for (int i=0; i < gsd->nbonds; ++i)
        {
        gsd->bond_from[i] = bonds[2*i] + 1;
        gsd->bond_to[i] = bonds[2*i+1] + 1;
        }
    SAFE_FREE(bonds);

    // successful read, so set pointers now
    *nbonds = gsd->nbonds;
    *from = gsd->bond_from;
    *to = gsd->bond_to;
    *nbondtypes = gsd->bondmap->ntypes;
    *bondtypename = gsd->bondmap->type;

    return MOLFILE_SUCCESS;
    }

//! Reads metadata about the GSD trajectory
/*!
 * \param mydata GSD trajectory
 * \param meta VMD trajectory metadata (output)
 *
 * The count of frames is set from the GSD trajectory. GSD almost always has
 * velocities stored in frame 0 unless the system was initialized to all 0s,
 * even if the user set them to "static" so that they aren't logged during the
 * simulation. Unfortunately, there's no way to tell if the intention is to use
 * the velocities from frame 0 (which may have been hacked to use as some coloring
 * field) statically, or if the intention is to omit them, and so velocities are
 * always read if we can find the chunk.
 *
 * \returns MOLFILE_SUCCESS
 */
static int read_gsd_timestep_metadata(void *mydata, molfile_timestep_metadata_t *meta)
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t *)mydata;

    meta->count = gsd->numframes;

    const struct gsd_index_entry* entry = gsd_find_chunk(gsd->handle, 0, "particles/velocity");
    meta->has_velocities = (entry != NULL);

    return MOLFILE_SUCCESS;
    }

//! Reads a timestep (frame) from the GSD trajectory
/*!
 * \param mydata GSD trajectory
 * \param natoms Number of atoms VMD expects in the frame
 * \param ts VMD timestep data (output)
 *
 * \returns MOLFILE_SUCCESS on success, MOLFILE_EOF when the last frame is read,
 *          and MOLFILE_ERROR on failure.
 *
 * If \a ts is a valid pointer and EOF has not been reached, the current frame is
 * read out of the GSD file. Because VMD only supports constant number of particles
 * (but GSD supports a changing number), the number of particles in the frame is
 * checked to ensure that it agrees with the number found in the first frame. An
 * error is reported if the number of particles changes.
 *
 * The timestep (integer) is reported as the "physical time" of the current frame.
 * The simulation box is also reported for each frame, and is converted from
 * HOOMD's triclinic tilt factors to the angles required by VMD. Particle positions
 * are read for each frame (an error is raised if the positions are not present).
 * The velocities are optionally (but in practice basically always) read from the
 * file, and default to frame 0 if necessary.
 *
 * If \a ts is NULL, then the frame is simply skipped. In both cases, the frame
 * counter is advanced in \a mydata.
 *
 * \note MOLFILE_EOF and MOLFILE_ERROR currently take the same value, so VMD
 *       cannot distinguish the two in the return value.
 */
static int read_gsd_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t*)mydata;

    if (gsd->frame >= gsd->numframes) return MOLFILE_EOF;

    if (ts != NULL)
        {
        // read the number of particles as a sanity check
        int cur_natoms = 0;
        int retval = read_chunk(gsd->handle, &cur_natoms, gsd->frame, "particles/N", sizeof(int), 0);
        if (retval == -3)
            {
            retval = read_chunk(gsd->handle, &cur_natoms, 0, "particles/N", sizeof(int), 0);
            }
        if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading number of particles from frame %d, aborting.\n", gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }
        else if (cur_natoms != natoms)
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) VMD does not support changing number of particles (%d in frame %d, but %d in frame 0), aborting.\n", cur_natoms, gsd->frame, natoms);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }

        // read frame timestep
        uint64_t timestep = 0;
        retval = read_chunk(gsd->handle, &timestep, gsd->frame, "configuration/step", sizeof(uint64_t), 0);
        if (retval == 0 || retval == -3)
            {
            ts->physical_time = (double)timestep;
            }
        else
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading timestep from frame %d, aborting.\n", gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }

        // read the box size, and convert tilt factors to angles
        float box[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f}; // default box specification
        retval = read_chunk(gsd->handle, box,  gsd->frame, "configuration/box", 6 * sizeof(float), 0);
        if (retval == -3) // extract from frame 0 otherwise
            {
            retval = read_chunk(gsd->handle, box, 0, "configuration/box", 6 * sizeof(float), 0);
            }
        // if retval is still nonzero, then there was an error, abort
        if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading box size from frame %d, aborting.\n", gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }
        else
            {
            ts->A = box[0]; ts->B = box[1]; ts->C = box[2];
            if (box[3] != 0.0f || box[4] != 0.0f || box[5] != 0.0f)
                {
                // need to resolve the tilt factors into angles
                const double xy = (double)box[3];
                const double xz = (double)box[4];
                const double yz = (double)box[5];
                const double norm1 = sqrt(1.0 + xy*xy);
                const double norm2 = sqrt(1.0 + xz*xz + yz*yz);

                const double cos_gamma= xy / norm1;
                const double cos_beta = xz / norm2;
                const double cos_alpha = (xy*xz + yz)/(norm1 * norm2);

                ts->alpha = (float)(acos(cos_alpha) * (180./M_PI));
                ts->beta = (float)(acos(cos_beta) * (180./M_PI));
                ts->gamma = (float)(acos(cos_gamma) * (180./M_PI));
                }
            else // orthorhombic
                {
                ts->alpha = 90.0f; ts->beta = 90.0f; ts->gamma = 90.0f;
                }
            }

        // read positions
        retval = read_chunk(gsd->handle, ts->coords, gsd->frame, "particles/position", 3*gsd->natoms*sizeof(float), gsd->natoms);
        if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading particle positions from frame %d, aborting.\n", gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }

        // read frame velocities
        if (ts->velocities != NULL)
            {
            retval = read_chunk(gsd->handle, ts->velocities, gsd->frame, "particles/velocity", 3*gsd->natoms*sizeof(float), gsd->natoms);
            if (retval == -3)
                {
                retval = read_chunk(gsd->handle, ts->velocities, 0, "particles/velocity", 3*gsd->natoms*sizeof(float), gsd->natoms);
                }

            if (retval != 0)
                {
                vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading particle velocities from frame %d, aborting.\n", gsd->frame);
                ++gsd->frame;
                return MOLFILE_ERROR;
                }
            }

        }
    ++gsd->frame;

    return MOLFILE_SUCCESS;
    }

//! Closes the GSD file for reading
/*!
 * \param mydata GSD trajectory
 *
 * \post All data stored in GSD trajectory is freed.
 *
 * \sa free_gsd_trajectory
 */
static void close_gsd_read(void *mydata)
    {
    free_gsd_trajectory((gsd_trajectory_t*)mydata);
    }

/* plugin registration */
//! Plugin object
static molfile_plugin_t plugin;

//! VMD plugin initialization
VMDPLUGIN_API int VMDPLUGIN_init()
    {
    memset(&plugin, 0, sizeof(molfile_plugin_t));
    plugin.abiversion = vmdplugin_ABIVERSION;
    plugin.type = MOLFILE_PLUGIN_TYPE;
    plugin.name = "gsd";
    plugin.prettyname = "HOOMD-blue GSD File";
    plugin.author = "Michael P. Howard";
    plugin.majorv = 0;
    plugin.minorv = 0;
    plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
    plugin.filename_extension = "gsd";

    plugin.open_file_read = open_gsd_read;
    plugin.read_structure = read_gsd_structure;
    plugin.read_bonds = read_gsd_bonds;
    plugin.read_next_timestep = read_gsd_timestep;
    plugin.read_timestep_metadata = read_gsd_timestep_metadata;
    plugin.close_file_read = close_gsd_read;
    return VMDPLUGIN_SUCCESS;
    }

//! VMD plugin registration
VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb)
    {
    (*cb)(v, (vmdplugin_t *)&plugin);
    return VMDPLUGIN_SUCCESS;
    }

//! VMD plugin finalization
VMDPLUGIN_API int VMDPLUGIN_fini()
    {
    return VMDPLUGIN_SUCCESS;
    }
