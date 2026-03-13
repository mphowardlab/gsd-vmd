// Copyright (c) 2017-2020, Michael P. Howard
// Copyright (c) 2021-2026, Auburn University
// Part of gsd-vmd, released under the BSD 3-Clause License.

#include "gsd.h"
#include "molfile_plugin.h"
#include "vmdconio.h"

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

static int safe_multiply_size(size_t N, size_t M, size_t element_size, size_t* total_size)
    {
    if (N == 0 || M == 0 || element_size == 0)
        {
        *total_size = 0;
        return 0;
        }

    if (N > SIZE_MAX / M)
        return -1;
    size_t num_elements = N * M;

    if (num_elements > SIZE_MAX / element_size)
        return -1;

    *total_size = num_elements * element_size;
    return 0;
    }

//! Safely allocate a chunk of 2D memory without allowing overflow
static void* safe_malloc(size_t N, size_t M, size_t element_size)
    {
    size_t num_bytes = 0;
    if (safe_multiply_size(N, M, element_size, &num_bytes) == 0 && num_bytes > 0)
        {
        return malloc(num_bytes);
        }
    else
        {
        return NULL;
        }
    }

//! Macro to safely free and NULL a pointer \a p
#define SAFE_FREE(p)  \
    do                \
        {             \
        if (p)        \
            {         \
            free(p);  \
            p = NULL; \
            }         \
        } while (0)

//! Resize memory
static int resize(void** buffer, size_t* buffer_capacity, size_t N, size_t M, size_t element_size)
    {
    size_t num_bytes = 0;
    if (safe_multiply_size(N, M, element_size, &num_bytes) != 0)
        return -1;

    // no need to resize, buffer always grows
    if (num_bytes <= *buffer_capacity)
        return 0;

    // reallocate the buffer
    SAFE_FREE(*buffer);
    *buffer = malloc(num_bytes);
    if (!*buffer)
        {
        *buffer_capacity = 0;
        return -1;
        }
    else
        {
        *buffer_capacity = num_bytes;
        return 0;
        }
    }

//! GSD handle object
typedef struct gsd_handle gsd_handle_t;

//! String type map
/*!
 * Maps an index to a string name
 */
typedef struct
    {
    int ntypes;  //!< Number of types mapped
    char** type; //!< Names of types
    } typemap_t;

//! Reallocate a type map and null entry strings
/*!
 * \param typemap String type map to reallocate
 * \param ntypes Number of types to allocate for
 * \returns 0 on success, -1 on failure
 *
 * When \a ntypes is 0 or the call fails, all memory is freed and pointers are nulled. When \a
 * ntypes is greater than 0, memory for the map is allocated and each entry initialized to a NULL
 * pointer when the call is successful.
 *
 * \warning Any values stored in the map are lost on reallocation.
 */
static int reallocate_typemap(typemap_t* typemap, int ntypes)
    {
    if (!typemap)
        return -1;

    // free any existing memory
    if (typemap->type)
        {
        for (int i = 0; i < typemap->ntypes; ++i)
            {
            SAFE_FREE(typemap->type[i]);
            }
        SAFE_FREE(typemap->type);
        }

    // try to allocate the memory for the requested number of types
    typemap->ntypes = ntypes;
    if (ntypes > 0)
        {
        typemap->type = (char**)safe_malloc(ntypes, 1, sizeof(char*));
        // quit early on failed allocation
        if (!typemap->type)
            {
            typemap->ntypes = 0;
            return -1;
            }

        // otherwise null every type pointer
        for (int i = 0; i < typemap->ntypes; ++i)
            {
            typemap->type[i] = NULL;
            }
        }

    return 0;
    }

//! Type map constructor
/*!
 * \param ntypes Number of types to hold in the map
 * \returns An initialized type map with all entries set to NULL pointers
 */
static typemap_t* allocate_typemap(int ntypes)
    {
    typemap_t* typemap = (typemap_t*)malloc(sizeof(typemap_t));
    if (!typemap)
        return NULL;

    // default initialize, then immediately reallocate to right size
    typemap->ntypes = 0;
    typemap->type = NULL;
    const int retval = reallocate_typemap(typemap, ntypes);
    if (retval != 0)
        {
        SAFE_FREE(typemap);
        }

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
    if (!typemap)
        return;

    // free any existing memory by calling a reallocation
    reallocate_typemap(typemap, 0);
    SAFE_FREE(typemap);
    }

//! GSD trajectory
typedef struct
    {
    gsd_handle_t* handle; //!< GSD file handle
    int frame;            //!< Current frame index
    int numframes;        //!< Number of frames in gsd file
    int natoms;           //!< Number of atoms in first frame
    typemap_t* typemap;   //!< Type map

    int nbonds;         //!< Number of bonds
    int* bond_from;     //!< First particle in bond (1-indexed)
    int* bond_to;       //!< Second particle in bond (2-indexed)
    typemap_t* bondmap; //!< Bond map
    } gsd_trajectory_t;

//! Constructor for GSD trajectory
/*!
 * \returns An allocated and default initialized GSD trajectory
 */
static gsd_trajectory_t* allocate_gsd_trajectory()
    {
    gsd_trajectory_t* gsd = (gsd_trajectory_t*)calloc(1, sizeof(gsd_trajectory_t));
    if (!gsd)
        return NULL;

    // try allocating memory first
    gsd->handle = (gsd_handle_t*)malloc(sizeof(gsd_handle_t));
    gsd->typemap = allocate_typemap(0);
    gsd->bondmap = allocate_typemap(0);
    if (!gsd->handle || !gsd->typemap || !gsd->bondmap)
        {
        SAFE_FREE(gsd->handle);
        free_typemap(gsd->typemap);
        free_typemap(gsd->bondmap);
        return NULL;
        }

    // finish default initializations
    gsd->frame = 0;
    gsd->numframes = 0;
    gsd->natoms = 0;

    gsd->nbonds = 0;
    gsd->bond_from = NULL;
    gsd->bond_to = NULL;

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
static void free_gsd_trajectory(gsd_trajectory_t* gsd)
    {
    if (!gsd)
        return;

    if (gsd->handle)
        {
        gsd_close(gsd->handle);
        }
    SAFE_FREE(gsd->handle);

    free_typemap(gsd->typemap);
    free_typemap(gsd->bondmap);

    SAFE_FREE(gsd->bond_from);
    SAFE_FREE(gsd->bond_to);

    SAFE_FREE(gsd);
    }

//! Read the size of a chunk element
static size_t read_chunk_element_size(gsd_handle_t* handle, uint64_t frame, const char* name)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, frame, name);
    return (entry) ? gsd_sizeof_type((enum gsd_type)entry->type) : 0;
    }

//! Read a chunk from the GSD file
/*!
 * \param handle Pointer to the GSD file handle
 * \param data Array to write the data into
 * \param frame Frame index
 * \param name Name of chunk to read
 * \param expected_N Expected number of elements in first dimension
 * \param expected_M Expected number of elements in second dimentions
 * \param element_size Expected bytes per element
 *
 * \returns gsd_error if entry exists, and 1 otherwise.
 *
 * If \a expected_N is nonzero, then the chunk size is validated to ensure it
 * contains \a expected_N entries.
 */
static int read_chunk(gsd_handle_t* handle,
                      void* data,
                      uint64_t frame,
                      const char* name,
                      size_t expected_N,
                      size_t expected_M,
                      size_t element_size)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, frame, name);
    if (!entry)
        {
        // silently ignore missing entries
        return 1;
        }
    else if (entry->N != expected_N || entry->M != expected_M
             || gsd_sizeof_type((enum gsd_type)entry->type) != element_size)
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) Incorrect shape of chunk '%s' at frame %d.\n",
                      name,
                      frame);
        return GSD_ERROR_FILE_CORRUPT;
        }

    int retval = gsd_read_chunk(handle, data, entry);
    if (retval == GSD_ERROR_IO)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Chunk '%s' : %s\n", name, strerror(errno));
        }
    else if (retval != GSD_SUCCESS)
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
static void* open_gsd_read(const char* filename, const char* filetype, int* natoms)
    {
    if (!filename)
        return NULL;

    gsd_trajectory_t* gsd = allocate_gsd_trajectory();
    if (!gsd)
        return NULL;

    int retval = gsd_open(gsd->handle, filename, GSD_OPEN_READONLY);
    if (retval == GSD_ERROR_IO)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) '%s': %s\n", filename, strerror(errno));
        }
    else if (retval == GSD_ERROR_NOT_A_GSD_FILE)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) '%s' is not a valid GSD file\n", filename);
        }
    else if (retval == GSD_ERROR_INVALID_GSD_FILE_VERSION)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Invalid GSD file version in '%s'\n", filename);
        }
    else if (retval == GSD_ERROR_FILE_CORRUPT)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Corrupt GSD file '%s'\n", filename);
        }
    else if (retval == GSD_ERROR_MEMORY_ALLOCATION_FAILED)
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) Unable to allocate memory opening '%s'\n",
                      filename);
        }
    else if (retval != GSD_SUCCESS)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error opening '%s'\n", filename);
        }
    // safe free and return in all cases of error
    if (retval != GSD_SUCCESS)
        {
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // validate schema
    if (strcmp(gsd->handle->header.schema, "hoomd") != 0)
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) Invalid schema in '%s', expecting 'hoomd'\n",
                      filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }
    if (gsd->handle->header.schema_version >= gsd_make_version(3, 0))
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) Invalid schema version in '%s', expecting >=1, <3\n",
                      filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // validate that at least one frame is written (is this just given?)
    gsd->numframes = gsd_get_nframes(gsd->handle);
    if (gsd->numframes == 0)
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) GSD file '%s' does not contain any frames\n",
                      filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // read the number of particles
    *natoms = 0;
    read_chunk(gsd->handle, natoms, 0, "particles/N", 1, 1, sizeof(int));
    if (*natoms == 0)
        {
        vmdcon_printf(VMDCON_ERROR,
                      "gsdplugin) No particles found in first frame of '%s'\n",
                      filename);
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
static int read_gsd_typemap(gsd_trajectory_t* gsd, molfile_atom_t* atoms)
    {
    if (!gsd || !atoms)
        return MOLFILE_ERROR;

    const struct gsd_index_entry* entry = gsd_find_chunk(gsd->handle, 0, "particles/types");
    if (entry) // types are present
        {
        char* data = (char*)safe_malloc(entry->N, entry->M, sizeof(char));
        if (!data)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error allocating buffer for particle type mapping\n");
            return MOLFILE_ERROR;
            }

        int retval = gsd_read_chunk(gsd->handle, data, entry);
        if (retval == GSD_ERROR_IO)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Type mapping 'particles/types' : %s\n",
                          strerror(errno));
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }
        else if (retval != GSD_SUCCESS)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error reading type mapping 'particles/types'\n");
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }

        // determine the maximum copy size from the molfile_atom_t
        const size_t max_name = sizeof(atoms->name);
        const size_t max_type = sizeof(atoms->name);
        const size_t max_nametype = (max_name < max_type) ? max_name - 1 : max_type - 1;
        if (max_nametype < (entry->M - 1))
            {
            vmdcon_printf(VMDCON_WARN,
                          "gsdplugin) Type names cannot exceed %d characters, truncating\n",
                          max_nametype);
            }

        // remalloc the type mapping and copy from the gsd data with null termination
        if (reallocate_typemap(gsd->typemap, entry->N) != 0)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error allocating memory for particle type mapping\n");
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }
        for (int i = 0; i < entry->N; ++i)
            {
            const char* name = data + i * entry->M;
            // get size of the name
            size_t l = strnlen(name, entry->M);
            if (l > max_nametype)
                {
                l = max_nametype;
                }

            gsd->typemap->type[i] = (char*)safe_malloc(l + 1, 1, sizeof(char));
            if (!gsd->typemap->type[i])
                {
                vmdcon_printf(VMDCON_ERROR,
                              "gsdplugin) Error allocating memory for particle type\n");
                reallocate_typemap(gsd->typemap, 0);
                SAFE_FREE(data);
                return MOLFILE_ERROR;
                }
            strncpy(gsd->typemap->type[i], name, l);
            gsd->typemap->type[i][l] = '\0'; // force null termination
            }
        SAFE_FREE(data);
        }
    else
        {
        // initialize default types from the HOOMD spec
        if (reallocate_typemap(gsd->typemap, 1) != 0)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error allocating memory for particle type mapping\n");
            return MOLFILE_ERROR;
            }
        gsd->typemap->type[0] = (char*)malloc(2 * sizeof(char));
        strncpy(gsd->typemap->type[0], "A\0", 2);
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
static int read_gsd_types(gsd_trajectory_t* gsd, molfile_atom_t* atoms)
    {
    if (!gsd || !atoms)
        return MOLFILE_ERROR;

    if (read_gsd_typemap(gsd, atoms) != MOLFILE_SUCCESS)
        return MOLFILE_ERROR;

    uint32_t* typeid = (uint32_t*)calloc(gsd->natoms, sizeof(uint32_t));
    if (!typeid)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error allocating buffer for partice typeids\n");
        return MOLFILE_ERROR;
        }

    int retval
        = read_chunk(gsd->handle, typeid, 0, "particles/typeid", gsd->natoms, 1, sizeof(uint32_t));
    if (retval == GSD_SUCCESS || retval == 1)
        {
        molfile_atom_t* a = atoms;

        for (int i = 0; i < gsd->natoms; ++i, ++a)
            {
            unsigned int typeid_i = typeid[i];

            if (typeid_i >= (unsigned int)gsd->typemap->ntypes)
                {
                vmdcon_printf(VMDCON_ERROR,
                              "gsdplugin) Invalid type ID %u for particle %d (max: %d)\n",
                              typeid_i,
                              i,
                              gsd->typemap->ntypes - 1);
                SAFE_FREE(typeid);
                return MOLFILE_ERROR;
                }

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
 * \param buffer Temporary memory
 * \param buffer_capacity Size of temporary memory
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int
read_gsd_mass(gsd_trajectory_t* gsd, molfile_atom_t* atoms, void** buffer, size_t* buffer_capacity)
    {
    if (!gsd || !atoms)
        return MOLFILE_ERROR;

    // resize buffer as needed
    const size_t element_size = read_chunk_element_size(gsd->handle, 0, "particles/mass");
    if (resize(buffer, buffer_capacity, gsd->natoms, 1, element_size) != 0)
        {
        return MOLFILE_ERROR;
        }

    // read values, accounting for different precisions
    int retval
        = read_chunk(gsd->handle, *buffer, 0, "particles/mass", gsd->natoms, 1, element_size);
    if (retval == GSD_SUCCESS)
        {
        molfile_atom_t* a = atoms;
        if (element_size == sizeof(double))
            {
            double* mass = (double*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->mass = mass[i];
                }
            }
        else
            {
            float* mass = (float*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->mass = mass[i];
                }
            }
        }
    else if (retval != 1)
        {
        return MOLFILE_ERROR;
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle charge from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 * \param buffer Temporary memory
 * \param buffer_capacity Size of temporary memory
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_charge(gsd_trajectory_t* gsd,
                           molfile_atom_t* atoms,
                           void** buffer,
                           size_t* buffer_capacity)
    {
    if (!gsd || !atoms)
        return MOLFILE_ERROR;

    // resize buffer as needed
    const size_t element_size = read_chunk_element_size(gsd->handle, 0, "particles/charge");
    if (resize(buffer, buffer_capacity, gsd->natoms, 1, element_size) != 0)
        {
        return MOLFILE_ERROR;
        }

    int retval
        = read_chunk(gsd->handle, *buffer, 0, "particles/charge", gsd->natoms, 1, element_size);
    if (retval == GSD_SUCCESS)
        {
        molfile_atom_t* a = atoms;
        if (element_size == sizeof(double))
            {
            double* charge = (double*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->charge = charge[i];
                }
            }
        else
            {
            float* charge = (float*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->charge = charge[i];
                }
            }
        }
    else if (retval != 1)
        {
        return MOLFILE_ERROR;
        }

    return MOLFILE_SUCCESS;
    }

//! Read particle radius from the GSD file
/*!
 * \param gsd GSD trajectory
 * \param atoms VMD atom properties
 * \param buffer Temporary memory
 * \param buffer_capacity Size of temporary memory
 *
 * \returns MOLFILE_SUCCESS on success or MOLFILE_ERROR on failure
 */
static int read_gsd_radius(gsd_trajectory_t* gsd,
                           molfile_atom_t* atoms,
                           void** buffer,
                           size_t* buffer_capacity)
    {
    if (!gsd || !atoms)
        return MOLFILE_ERROR;

    // resize buffer as needed
    const size_t element_size = read_chunk_element_size(gsd->handle, 0, "particles/diameter");
    if (resize(buffer, buffer_capacity, gsd->natoms, 1, element_size) != 0)
        {
        return MOLFILE_ERROR;
        }

    int retval
        = read_chunk(gsd->handle, *buffer, 0, "particles/diameter", gsd->natoms, 1, element_size);
    if (retval == GSD_SUCCESS)
        {
        molfile_atom_t* a = atoms;
        if (element_size == sizeof(double))
            {
            double* diameter = (double*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->radius = 0.5f * diameter[i];
                }
            }
        else
            {
            float* diameter = (float*)*buffer;
            for (int i = 0; i < gsd->natoms; ++i, ++a)
                {
                a->radius = 0.5f * diameter[i];
                }
            }
        }
    else if (retval != 1)
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
static int read_gsd_structure(void* mydata, int* optflags, molfile_atom_t* atoms)
    {
    gsd_trajectory_t* gsd = (gsd_trajectory_t*)mydata;

        // loop through the atoms and fill in with defaults per the HOOMD GSD schema
        {
        molfile_atom_t* a = atoms;
        for (int i = 0; i < gsd->natoms; ++i, ++a)
            {
            strncpy(a->name, "A\0", sizeof(atoms->name));
            strncpy(a->type, "A\0", sizeof(atoms->type));
            a->resname[0] = '\0';
            a->resid = 0;
            a->segid[0] = '\0';
            a->chain[0] = '\0';

            a->mass = 1.0f;
            a->charge = 0.0f;
            a->radius = 0.5f;
            }
        }
    *optflags = MOLFILE_MASS | MOLFILE_CHARGE | MOLFILE_RADIUS;

    // map the particle types
    if (read_gsd_types(gsd, atoms) != MOLFILE_SUCCESS)
        return MOLFILE_ERROR;

    // shared buffer for optional properties
    void* buffer = NULL;
    size_t buffer_capacity = 0;

    // mass
    if (read_gsd_mass(gsd, atoms, &buffer, &buffer_capacity) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(buffer);
        return MOLFILE_ERROR;
        }
    // charge
    if (read_gsd_charge(gsd, atoms, &buffer, &buffer_capacity) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(buffer);
        return MOLFILE_ERROR;
        }
    // radius
    if (read_gsd_radius(gsd, atoms, &buffer, &buffer_capacity) != MOLFILE_SUCCESS)
        {
        SAFE_FREE(buffer);
        return MOLFILE_ERROR;
        }
    SAFE_FREE(buffer);

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
static int read_bondmap(gsd_handle_t* handle, const char* name, typemap_t* bondmap)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, 0, name);
    if (entry) // types are present
        {
        char* data = (char*)safe_malloc(entry->N, entry->M, sizeof(char));
        int retval = gsd_read_chunk(handle, data, entry);
        if (retval == GSD_ERROR_IO)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Bond type mapping '%s' : %s\n",
                          name,
                          strerror(errno));
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }
        else if (retval != GSD_SUCCESS)
            {
            vmdcon_printf(VMDCON_ERROR, "gsdplugin) Error reading bond type mapping '%s'\n", name);
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }

        if (reallocate_typemap(bondmap, entry->N) != 0)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error allocating memory for bond type mapping\n");
            SAFE_FREE(data);
            return MOLFILE_ERROR;
            }
        if (bondmap && bondmap->type)
            {
            for (int i = 0; i < entry->N; ++i)
                {
                const char* name = data + i * entry->M;
                // get size of the name
                size_t l = strnlen(name, entry->M);

                // resizing guarantees that all member chars are nulled, so can just malloc
                bondmap->type[i] = (char*)malloc((l + 1) * sizeof(char));
                if (!bondmap->type[i])
                    {
                    vmdcon_printf(VMDCON_ERROR,
                                  "gsdplugin) Error allocating memory for particle type\n");
                    reallocate_typemap(bondmap, 0);
                    SAFE_FREE(data);
                    return MOLFILE_ERROR;
                    }
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
static int read_gsd_bonds(void* mydata,
                          int* nbonds,
                          int** from,
                          int** to,
                          float** bondorder,
                          int** bondtype,
                          int* nbondtypes,
                          char*** bondtypename)
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
    gsd_trajectory_t* gsd = (gsd_trajectory_t*)mydata;
    int retval = read_chunk(gsd->handle, &gsd->nbonds, 0, "bonds/N", 1, 1, sizeof(int));
    if (retval == 1 || gsd->nbonds == 0)
        {
        // return successfully if bonds are not present in the file
        return MOLFILE_SUCCESS;
        }
    else if (retval != GSD_SUCCESS)
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
    uint32_t* bonds = (uint32_t*)safe_malloc(gsd->nbonds, 2, sizeof(uint32_t));
    gsd->bond_from = (int*)safe_malloc(gsd->nbonds, 1, sizeof(int));
    gsd->bond_to = (int*)safe_malloc(gsd->nbonds, 1, sizeof(int));
    if (!bonds || !gsd->bond_from || !gsd->bond_to)
        {
        vmdcon_printf(VMDCON_ERROR, "gsdplugin) Failed to allocate memory for bonds\n");
        SAFE_FREE(bonds);
        SAFE_FREE(gsd->bond_from);
        SAFE_FREE(gsd->bond_to);
        return MOLFILE_ERROR;
        }

    retval = read_chunk(gsd->handle, bonds, 0, "bonds/group", gsd->nbonds, 2, sizeof(uint32_t));
    if (retval != GSD_SUCCESS)
        {
        SAFE_FREE(bonds);
        SAFE_FREE(gsd->bond_from);
        SAFE_FREE(gsd->bond_to);
        return MOLFILE_ERROR;
        }
    for (int i = 0; i < gsd->nbonds; ++i)
        {
        gsd->bond_from[i] = bonds[2 * i] + 1;
        gsd->bond_to[i] = bonds[2 * i + 1] + 1;
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
static int read_gsd_timestep_metadata(void* mydata, molfile_timestep_metadata_t* meta)
    {
    gsd_trajectory_t* gsd = (gsd_trajectory_t*)mydata;

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
static int read_gsd_timestep(void* mydata, int natoms, molfile_timestep_t* ts)
    {
    gsd_trajectory_t* gsd = (gsd_trajectory_t*)mydata;

    if (gsd->frame >= gsd->numframes)
        return MOLFILE_EOF;

    if (ts)
        {
        // read the number of particles as a sanity check
        int cur_natoms = 0;
        int retval
            = read_chunk(gsd->handle, &cur_natoms, gsd->frame, "particles/N", 1, 1, sizeof(int));
        if (retval == 1)
            {
            retval = read_chunk(gsd->handle, &cur_natoms, 0, "particles/N", 1, 1, sizeof(int));
            }
        if (retval != GSD_SUCCESS)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error reading number of particles from frame %d, aborting.\n",
                          gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }
        else if (cur_natoms != natoms)
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) VMD does not support changing number of particles (%d in "
                          "frame %d, but %d in frame 0), aborting.\n",
                          cur_natoms,
                          gsd->frame,
                          natoms);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }

        // read frame timestep
        uint64_t timestep = 0;
        retval = read_chunk(gsd->handle,
                            &timestep,
                            gsd->frame,
                            "configuration/step",
                            1,
                            1,
                            sizeof(uint64_t));
        if (retval == GSD_SUCCESS || retval == 1)
            {
            ts->physical_time = (double)timestep;
            }
        else
            {
            vmdcon_printf(VMDCON_ERROR,
                          "gsdplugin) Error reading timestep from frame %d, aborting.\n",
                          gsd->frame);
            ++gsd->frame;
            return MOLFILE_ERROR;
            }

        // read the box and convert tilt factors to angles
        float box[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f}; // default box specification
            {
            // get size of float for box, using either current frame or frame 0
            uint64_t box_frame = gsd->frame;
            size_t element_size
                = read_chunk_element_size(gsd->handle, box_frame, "configuration/box");
            if (element_size == 0)
                {
                box_frame = 0;
                read_chunk_element_size(gsd->handle, box_frame, "configuration/box");
                }

            int retval = GSD_ERROR_FILE_CORRUPT;
            if (element_size == sizeof(double))
                {
                double buffer[6];
                retval = read_chunk(gsd->handle,
                                    buffer,
                                    box_frame,
                                    "configuration/box",
                                    6,
                                    1,
                                    element_size);
                if (retval == GSD_SUCCESS)
                    {
                    for (int i = 0; i < 6; ++i)
                        {
                        box[i] = buffer[i];
                        }
                    }
                }
            else if (element_size == sizeof(float))
                {
                retval = read_chunk(gsd->handle,
                                    box,
                                    box_frame,
                                    "configuration/box",
                                    6,
                                    1,
                                    element_size);
                }

            if (retval != GSD_SUCCESS)
                {
                vmdcon_printf(VMDCON_ERROR,
                              "gsdplugin) Error reading box from frame %d, aborting.\n",
                              gsd->frame);
                ++gsd->frame;
                return MOLFILE_ERROR;
                }

            if (box[3] != 0.0f || box[4] != 0.0f || box[5] != 0.0f)
                {
                // define lattice constants in terms of box size and tilt factors
                const double xy = (double)box[3];
                const double xz = (double)box[4];
                const double yz = (double)box[5];
                const double norm1 = sqrt(1.0 + xy * xy);
                const double norm2 = sqrt(1.0 + xz * xz + yz * yz);

                ts->A = box[0];
                ts->B = (float)(norm1 * box[1]);
                ts->C = (float)(norm2 * box[2]);

                // need to resolve the tilt factors into angles
                const double cos_gamma = xy / norm1;
                const double cos_beta = xz / norm2;
                const double cos_alpha = (xy * xz + yz) / (norm1 * norm2);

                ts->alpha = (float)(acos(cos_alpha) * (180. / M_PI));
                ts->beta = (float)(acos(cos_beta) * (180. / M_PI));
                ts->gamma = (float)(acos(cos_gamma) * (180. / M_PI));
                }
            else // orthorhombic
                {
                ts->A = box[0];
                ts->B = box[1];
                ts->C = box[2];
                ts->alpha = 90.0f;
                ts->beta = 90.0f;
                ts->gamma = 90.0f;
                }
            }

        // buffer for positions and velocities
        double* particle_buffer = NULL;
        size_t particle_buffer_capacity = 0;

            // read positions
            {
            size_t element_size
                = read_chunk_element_size(gsd->handle, gsd->frame, "particles/position");

            int retval = GSD_ERROR_FILE_CORRUPT;
            if (element_size == sizeof(double))
                {
                if (resize((void**)&particle_buffer,
                           &particle_buffer_capacity,
                           gsd->natoms,
                           3,
                           element_size)
                    != 0)
                    {
                    vmdcon_printf(VMDCON_ERROR,
                                  "gsdplugin) Error allocating buffer for particle positions from "
                                  "frame %d, aborting.\n",
                                  gsd->frame);
                    SAFE_FREE(particle_buffer);
                    ++gsd->frame;
                    return MOLFILE_ERROR;
                    }

                retval = read_chunk(gsd->handle,
                                    particle_buffer,
                                    gsd->frame,
                                    "particles/position",
                                    gsd->natoms,
                                    3,
                                    element_size);

                if (retval == GSD_SUCCESS)
                    {
                    // this call is always safe because the malloc above succeed
                    size_t num_elements = 0;
                    safe_multiply_size(gsd->natoms, 3, 1, &num_elements);
                    for (size_t i = 0; i < num_elements; ++i)
                        {
                        ts->coords[i] = particle_buffer[i];
                        }
                    }
                }
            else if (element_size == sizeof(float))
                {
                retval = read_chunk(gsd->handle,
                                    ts->coords,
                                    gsd->frame,
                                    "particles/position",
                                    gsd->natoms,
                                    3,
                                    element_size);
                }

            if (retval != GSD_SUCCESS)
                {
                vmdcon_printf(
                    VMDCON_ERROR,
                    "gsdplugin) Error reading particle positions from frame %d, aborting.\n",
                    gsd->frame);
                SAFE_FREE(particle_buffer);
                ++gsd->frame;
                return MOLFILE_ERROR;
                }
            }

        // read frame velocities
        if (ts->velocities != NULL)
            {
            uint64_t velocity_frame = gsd->frame;
            size_t element_size
                = read_chunk_element_size(gsd->handle, velocity_frame, "particles/velocity");
            if (element_size == 0)
                {
                velocity_frame = 0;
                element_size
                    = read_chunk_element_size(gsd->handle, velocity_frame, "particles/velocity");
                }

            int retval = GSD_ERROR_FILE_CORRUPT;
            if (element_size == sizeof(double))
                {
                if (resize((void**)&particle_buffer,
                           &particle_buffer_capacity,
                           gsd->natoms,
                           3,
                           element_size)
                    != 0)
                    {
                    vmdcon_printf(VMDCON_ERROR,
                                  "gsdplugin) Error allocating buffer for particle velocities from "
                                  "frame %d, aborting.\n",
                                  gsd->frame);
                    SAFE_FREE(particle_buffer);
                    ++gsd->frame;
                    return MOLFILE_ERROR;
                    }

                retval = read_chunk(gsd->handle,
                                    particle_buffer,
                                    velocity_frame,
                                    "particles/velocity",
                                    gsd->natoms,
                                    3,
                                    element_size);

                if (retval == GSD_SUCCESS)
                    {
                    // this call is always safe because the malloc above succeed
                    size_t num_elements = 0;
                    safe_multiply_size(gsd->natoms, 3, 1, &num_elements);
                    for (size_t i = 0; i < num_elements; ++i)
                        {
                        ts->velocities[i] = particle_buffer[i];
                        }
                    }
                }
            else if (element_size == sizeof(float))
                {
                retval = read_chunk(gsd->handle,
                                    ts->velocities,
                                    velocity_frame,
                                    "particles/velocity",
                                    gsd->natoms,
                                    3,
                                    element_size);
                }

            if (retval != GSD_SUCCESS)
                {
                vmdcon_printf(
                    VMDCON_ERROR,
                    "gsdplugin) Error reading particle velocities from frame %d, aborting.\n",
                    gsd->frame);
                SAFE_FREE(particle_buffer);
                ++gsd->frame;
                return MOLFILE_ERROR;
                }
            }

        SAFE_FREE(particle_buffer);
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
static void close_gsd_read(void* mydata)
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
    plugin.minorv = 5;
#ifdef _WIN32
    plugin.is_reentrant = VMDPLUGIN_THREADUNSAFE;
#else
    plugin.is_reentrant = VMDPLUGIN_THREADSAFE;
#endif
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
VMDPLUGIN_API int VMDPLUGIN_register(void* v, vmdplugin_register_cb cb)
    {
    (*cb)(v, (vmdplugin_t*)&plugin);
    return VMDPLUGIN_SUCCESS;
    }

//! VMD plugin finalization
VMDPLUGIN_API int VMDPLUGIN_fini()
    {
    return VMDPLUGIN_SUCCESS;
    }
