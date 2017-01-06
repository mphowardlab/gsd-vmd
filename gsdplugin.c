#include "gsd.h"
#include "molfile_plugin.h"
#include "vmdconio.h"

#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

//! GSD trajectory
typedef struct
    {
    struct gsd_handle handle;   //!< GSD file handle
    int frame;                  //!< Current frame index
    int numframes;              //!< Number of frames in gsd file
    int natoms;                 //!< Number of atoms in first frame
    int numtypes;               //!< Number of types
    char **typemap;             //!< Names of types

    int nbonds;                 //!< Number of bonds
    int *bond_from;             //!< First particle in bond (1-indexed)
    int *bond_to;               //!< Second particle in bond (2-indexed)
    int numbondtypes;           //!< Number of bond types
    char **bondmap;             //!< Names of bonds
    } gsd_trajectory_t;

static void free_gsd_typemap(gsd_trajectory_t *gsd)
    {
    // free the typemap memory if it has been allocated
    for (int i=0; i < gsd->numtypes; ++i)
        {
        if (gsd->typemap[i]) free(gsd->typemap[i]);
        gsd->typemap[i] = NULL;
        }
    if (gsd->typemap) free(gsd->typemap);
    gsd->typemap = NULL;
    }

static void free_gsd_bonds(gsd_trajectory_t *gsd)
    {
    gsd->nbonds = 0;
    if (gsd->bond_from)
        {
        free(gsd->bond_from);
        gsd->bond_from = NULL;
        }
    if (gsd->bond_to)
        {
        free(gsd->bond_to);
        gsd->bond_to = NULL;
        }
    }

static void free_gsd_bondmap(gsd_trajectory_t *gsd)
    {
    gsd->numbondtypes = 0;
    if (gsd->bondmap)
        {
        free(gsd->bondmap);
        gsd->bondmap = NULL;
        }
    }

//! Destructor for gsd_trajectory_t
static void free_gsd_trajectory(gsd_trajectory_t *gsd)
    {
    if (gsd)
        {
        gsd_close(&gsd->handle);
        free_gsd_typemap(gsd);
        free_gsd_bonds(gsd);
        free_gsd_bondmap(gsd);
        free(gsd);
        }
    gsd = NULL;
    }

//! Read a chunk from the GSD file
/*!
 * \param handle Pointer to the GSD file handle
 * \param data Array to write the data into
 * \param frame Frame index
 * \param name Name of chunk to read
 * \param expected_size Expected number of elements in the chunk
 * \param expected_N Expected number of particles
 */
static int read_chunk(struct gsd_handle *handle,
                       void *data,
                       uint64_t frame,
                       const char *name,
                       size_t expected_size,
                       unsigned int expected_N,
                       bool log_console)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(handle, frame, name);
    if (entry == NULL)
        {
        // silently ignore missing entries
        return -3;
        }
    else if (expected_N != 0 && entry->N != expected_N)
        {
        if (log_console)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Incorrect number of entries in chunk '%s' at frame %d.\n"
                                        "     Expected %d and found %d\n", name, frame, expected_N, entry->N);
            }
        return -2;
        }

    // check that size of entry matches the expected size
    size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
    if (actual_size != expected_size)
        {
        if (log_console)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Incorrect entry size in chunk '%s' at frame %d.\n"
                                        "     Expected %d B and found %d B\n", name, frame, expected_size, actual_size);
            }
        return false;
        }

    int retval = gsd_read_chunk(handle, data, entry);
    if (log_console)
        {
        if (retval == -1)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Chunk '%s' : %s\n", name, strerror(errno));
            }
        else if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Error reading chunk '%s'\n", name);
            }
        }

    return retval;
    }

static void *open_gsd_read(const char *filename, const char *filetype, int *natoms)
    {
    if (!filename) return NULL;

    gsd_trajectory_t *gsd = (gsd_trajectory_t *)calloc(1, sizeof(gsd_trajectory_t));
    if (!gsd) return NULL;

    int retval = gsd_open(&gsd->handle, filename, GSD_OPEN_READONLY);
    if (retval == -1)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) '%s': %s\n", filename, strerror(errno));
        }
    else if (retval == -2)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) '%s' is not a valid GSD file\n", filename);
        }
    else if (retval == -3)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) Invalid GSD file version in '%s'\n", filename);
        }
    else if (retval == -4)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) Corrupt GSD file '%s'\n", filename);
        }
    else if (retval == -5)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) Unable to allocate memory opening '%s'\n", filename);
        }
    else if (retval != 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) Error opening '%s'\n", filename);
        }
    // safe free and return in all cases of error
    if (retval != 0)
        {
        free_gsd_trajectory(gsd);
        return NULL;
        }

    // validate schema

    // validate that at least one frame is written (is this just given?)
    gsd->numframes = gsd_get_nframes(&gsd->handle);

    // read the number of particles
    *natoms = 0;
    read_chunk(&gsd->handle, natoms, 0, "particles/N", 4, 0, false);
    if (*natoms == 0)
        {
        vmdcon_printf(VMDCON_ERROR, "gsd) No particles found in first frame of '%s'\n", filename);
        free_gsd_trajectory(gsd);
        return NULL;
        }

    gsd->natoms = *natoms;
    gsd->frame = 0;
    gsd->numtypes = 0;
    gsd->typemap = NULL;

    gsd->nbonds = 0;
    gsd->bond_to = NULL;
    gsd->bond_from = NULL;
    gsd->numbondtypes = 0;
    gsd->bondmap = NULL;

    return gsd;
    }

static int read_typemap(gsd_trajectory_t *gsd, molfile_atom_t *atoms)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(&gsd->handle, 0, "particles/types");
    if (entry != NULL) // types are present
        {
        size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
        char* data = (char*)malloc(actual_size);
        int retval = gsd_read_chunk(&gsd->handle, data, entry);
        if (retval == -1)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Type mapping 'particles/types' : %s\n", strerror(errno));
            return MOLFILE_ERROR;
            }
        else if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Error reading type mapping 'particles/types'\n");
            return MOLFILE_ERROR;
            }

        // determine the maximum copy size from the molfile_atom_t
        if (atoms == NULL) return MOLFILE_ERROR;
        const size_t max_name = sizeof(atoms->name);
        const size_t max_type = sizeof(atoms->name);
        const size_t max_nametype = (max_name < max_type) ? max_name-1 : max_type-1;
        if (max_nametype < (entry->M - 1))
            {
            vmdcon_printf(VMDCON_WARN, "gsd) Type names cannot exceed %d characters, truncating\n", max_nametype);
            }

        // remalloc the type mapping and copy from the gsd data with null termination
        free_gsd_typemap(gsd);
        gsd->numtypes = entry->N;
        gsd->typemap = (char**)malloc(gsd->numtypes * sizeof(char*));
        for (int i=0; i < gsd->numtypes; ++i)
            {
            const char *name = data + i*entry->M;
            // get size of the name
            size_t l = strnlen(name, entry->M);
            if (l > max_nametype)
                {
                l = max_nametype;
                }

            gsd->typemap[i] = (char*)malloc((l+1) * sizeof(char));
            strncpy(gsd->typemap[i], name, l);
            gsd->typemap[i][l] = '\0'; // force null termination
            }
        }
    else
        {
        // initialize default types from the HOOMD spec
        gsd->numtypes = 1;
        gsd->typemap = (char**)malloc(sizeof(char*));
        gsd->typemap[0] = (char*)malloc(2*sizeof(char));
        strncpy(gsd->typemap[0],"A\0",2);
        }

    return MOLFILE_SUCCESS;
    }

static int read_bondmap(gsd_trajectory_t *gsd, const char *name, int *numbondtypes, char ***bondmap)
    {
    const struct gsd_index_entry* entry = gsd_find_chunk(&gsd->handle, 0, name);
    if (entry != NULL) // types are present
        {
        size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);
        char* data = (char*)malloc(actual_size);
        int retval = gsd_read_chunk(&gsd->handle, data, entry);
        if (retval == -1)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Type mapping '%s' : %s\n", name, strerror(errno));
            return MOLFILE_ERROR;
            }
        else if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Error reading type mapping '%s'\n", name);
            return MOLFILE_ERROR;
            }

        *numbondtypes = entry->N;
        *bondmap = (char**)malloc(gsd->numtypes * sizeof(char*));
        for (int i=0; i < entry->N; ++i)
            {
            const char *name = data + i*entry->M;
            // get size of the name
            size_t l = strnlen(name, entry->M);

            (*bondmap)[i] = (char*)malloc((l+1) * sizeof(char));
            strncpy((*bondmap)[i], name, l);
            (*bondmap)[i][l] = '\0'; // force null termination
            }
        }

    return MOLFILE_SUCCESS;
    }

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
        {
        if (read_typemap(gsd, atoms) != MOLFILE_SUCCESS) return MOLFILE_ERROR;
        uint32_t *typeid = (uint32_t*)calloc(gsd->natoms, sizeof(uint32_t));
        int retval = read_chunk(&gsd->handle, typeid, 0, "particles/typeid", gsd->natoms * sizeof(uint32_t), gsd->natoms, true);
        if (retval == 0 || retval == -3)
            {
            molfile_atom_t *a = atoms;
            for (int i=0; i < gsd->natoms; ++i, ++a)
                {
                unsigned int typeid_i = typeid[i];
                strncpy(a->name, gsd->typemap[typeid_i], sizeof(atoms->name));
                strncpy(a->type, gsd->typemap[typeid_i], sizeof(atoms->type));
                }
            }
        else
            {
            if (typeid)
                {
                free(typeid);
                typeid = NULL;
                }
            return MOLFILE_ERROR;
            }
        if (typeid)
            {
            free(typeid);
            typeid = NULL;
            }
        }

    // try to set optional properties
    float *props = (float*)calloc(gsd->natoms, sizeof(float));
    // mass
        {
        int retval = read_chunk(&gsd->handle, props, 0, "particles/mass", gsd->natoms * sizeof(float), gsd->natoms, true);
        if (retval == 0)
            {
            molfile_atom_t *a = atoms;
            for (int i=0; i < gsd->natoms; ++i, ++a)
                {
                a->mass = props[i];
                }
            }
        else if (retval != -3)
            {
            if (props)
                {
                free(props);
                props = NULL;
                }
            return MOLFILE_ERROR;
            }
        }
    // charge
        {
        int retval = read_chunk(&gsd->handle, props, 0, "particles/charge", gsd->natoms * sizeof(float), gsd->natoms, true);
        if (retval == 0)
            {
            molfile_atom_t *a = atoms;
            for (int i=0; i < gsd->natoms; ++i, ++a)
                {
                a->charge = props[i];
                }
            }
        else if (retval != -3)
            {
            if (props)
                {
                free(props);
                props = NULL;
                }
            return MOLFILE_ERROR;
            }
        }
    // radius
        {
        int retval = read_chunk(&gsd->handle, props, 0, "particles/diameter", gsd->natoms * sizeof(float), gsd->natoms, true);
        if (retval == 0)
            {
            molfile_atom_t *a = atoms;
            for (int i=0; i < gsd->natoms; ++i, ++a)
                {
                a->radius = 0.5f * props[i];
                }
            }
        else if (retval != -3)
            {
            if (props)
                {
                free(props);
                props = NULL;
                }
            return MOLFILE_ERROR;
            }
        }
    if (props)
        {
        free(props);
        props = NULL;
        }

    return MOLFILE_SUCCESS;
    }

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

    gsd_trajectory_t *gsd = (gsd_trajectory_t*)mydata;

    // default is to not supply any bonds
    *nbonds = 0;
    *from = NULL;
    *to = NULL;
    *bondtype = NULL;
    *nbondtypes = 0;
    *bondtypename = NULL;

    // check number of bonds and exit early if no bonds are present, or on read error
    int retval = read_chunk(&gsd->handle, &gsd->nbonds,  0, "bonds/N", sizeof(int), 0, true);
    if (retval == -3 || gsd->nbonds == 0)
        {
        // return without reading bonds if are not present
        return MOLFILE_SUCCESS;
        }
    else if (retval != 0)
        {
        // this is a read error, quit early
        return MOLFILE_ERROR;
        }

    // acquire the bondname map
    if (read_bondmap(gsd, "bonds/types", &gsd->numbondtypes, &gsd->bondmap) != MOLFILE_SUCCESS)
        {
        free_gsd_bondmap(gsd);
        return MOLFILE_ERROR;
        }

    // read in the bonds, and remap them with 1-indexing
    uint32_t *bonds = (uint32_t*)calloc(2*gsd->nbonds, sizeof(uint32_t));
    retval = read_chunk(&gsd->handle, bonds, 0, "bonds/group", 2*gsd->nbonds*sizeof(uint32_t), gsd->nbonds, true);
    if (retval != 0)
        {
        free(bonds);
        bonds = NULL;
        return MOLFILE_ERROR;
        }
    gsd->bond_from = (int*)malloc(gsd->nbonds * sizeof(int));
    gsd->bond_to = (int*)malloc(gsd->nbonds * sizeof(int));
    for (int i=0; i < gsd->nbonds; ++i)
        {
        gsd->bond_from[i] = bonds[2*i] + 1;
        gsd->bond_to[i] = bonds[2*i+1] + 1;
        }
    // safe-free of the temporary data now that we are done with it
    if (bonds)
        {
        free(bonds);
        bonds = NULL;
        }

    // successful read, so set pointers now
    *nbonds = gsd->nbonds;
    *from = gsd->bond_from;
    *to = gsd->bond_to;
    *nbondtypes = gsd->numbondtypes;
    *bondtypename = gsd->bondmap;

    return MOLFILE_SUCCESS;
    }

static int read_gsd_timestep_metadata(void *mydata, molfile_timestep_metadata_t *meta)
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t *)mydata;

    meta->count = gsd->numframes;

    const struct gsd_index_entry* entry = gsd_find_chunk(&gsd->handle, 0, "particles/velocity");
    meta->has_velocities = (entry != NULL);

    return MOLFILE_SUCCESS;
    }

static int read_gsd_timestep(void *mydata, int natoms, molfile_timestep_t *ts)
    {
    gsd_trajectory_t *gsd = (gsd_trajectory_t*)mydata;

    if (gsd->frame >= gsd->numframes) return MOLFILE_EOF;

    if (ts != NULL)
        {
        // read frame timestep
        uint64_t timestep = 0;
        retval = read_chunk(&gsd->handle, &timestep, gsd->frame, "configuration/step", sizeof(uint64_t), 0, true);
        if (retval == 0 || retval == -3)
            {
            ts->physical_time = (double)timestep;
            }
        else
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Error reading timestep from frame %d', aborting.\n", gsd->frame);
            return MOLFILE_ERROR;
            }

        float box[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f}; // default box specification
        int retval = read_chunk(&gsd->handle, box,  gsd->frame, "configuration/box", 6 * sizeof(float), 0, true);
        if (retval == -3) // extract from frame 0 otherwise
            {
            retval = read_chunk(&gsd->handle, box, 0, "configuration/box", 6 * sizeof(float), 0, true);
            }
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

        // read positions
        retval = read_chunk(&gsd->handle, ts->coords, gsd->frame, "particles/position", 3*gsd->natoms*sizeof(float), gsd->natoms, true);
        if (retval != 0)
            {
            vmdcon_printf(VMDCON_ERROR, "gsd) Error reading particle positions from frame %d', aborting.\n", gsd->frame);
            return MOLFILE_ERROR;
            }

        // read frame velocities
        if (ts->velocities != NULL)
            {
            retval = read_chunk(&gsd->handle, ts->velocities, gsd->frame, "particles/velocity", 3*gsd->natoms*sizeof(float), gsd->natoms, true);
            if (retval == -3)
                {
                retval = read_chunk(&gsd->handle, ts->velocities, 0, "particles/velocity", 3*gsd->natoms*sizeof(float), gsd->natoms, true);
                }
            
            if (retval != 0)
                {
                vmdcon_printf(VMDCON_ERROR, "gsd) Error reading particle velocities from frame %d, aborting.\n", gsd->frame);
                return MOLFILE_ERROR;
                }
            }

        }
    ++gsd->frame;

    return MOLFILE_SUCCESS;
    }

static void close_gsd_read(void *mydata)
    {
    free_gsd_trajectory(mydata);
    }

/* plugin registration */
static molfile_plugin_t plugin;

VMDPLUGIN_API int VMDPLUGIN_init() {
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

VMDPLUGIN_API int VMDPLUGIN_register(void *v, vmdplugin_register_cb cb)
    {
    (*cb)(v, (vmdplugin_t *)&plugin);
    return VMDPLUGIN_SUCCESS;
    }

VMDPLUGIN_API int VMDPLUGIN_fini()
    {
    return VMDPLUGIN_SUCCESS;
    }
