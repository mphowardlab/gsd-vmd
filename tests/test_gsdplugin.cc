#include "catch.hpp"
#include "molfile_plugin.h"

#include <dlfcn.h>
#include <iostream>
#include <string>

extern "C"
    {
    typedef int (*initfunc)(void);
    typedef int (*regfunc)(void *, vmdplugin_register_cb);
    typedef int (*finifunc)(void);
    }

//! Opens a dynamic library on linux
class DynamicLibrary
    {
    public:
        DynamicLibrary()
            : handle_(nullptr)
            {
            }

        //! Constructor
        DynamicLibrary(const std::string& lib)
            : handle_(nullptr), lib_(lib)
            {
            }

        //! Destructor
        ~DynamicLibrary()
            {
            close();
            }

        //! Opens the library
        /*!
         * \return True if library is opened, false on failure
         */
        bool open()
            {
            // don't open twice
            if (handle_) return true;

            // flush the error handler
            dlerror();
            handle_ = dlopen(lib_.c_str(), RTLD_NOW);
            if (!handle_)
                {
                std::cerr << dlerror() << std::endl;
                handle_ = nullptr;
                return false;
                }
            else
                {
                return true;
                }
            }

        //! Load a symbol from the library
        void* load(const std::string& name) const
            {
            if (!handle_) return nullptr;

            // flush error handling
            dlerror();

            void *result = dlsym(handle_, name.c_str());
            if (!result)
                {
                std::cerr << dlerror() << std::endl;
                return nullptr;
                }
            return result;
            }

        //! Close the library
        void close()
            {
            if (handle_)
                {
                dlclose(handle_);
                handle_ = nullptr;
                }
            }

    private:
        void *handle_; //!< Handle to the library
        std::string lib_;   //!< Name of the library
    };

//! GSD plugin wrapper
class GSDPlugin
    {
    public:
        //! Constructor
        GSDPlugin()
            : plugin_(nullptr), init_(nullptr), fini_(nullptr)
            {
            lib = DynamicLibrary("./gsdplugin.so");

            // initialize plugin, and call register to get a pointer to the plugin
            if (init())
                {
                void *rfunc = lib.load("vmdplugin_register");
                if (rfunc)
                    {
                    ((regfunc)rfunc)(this, setPlugin);
                    }
                }
            }

        //! Accessor for the plugin
        molfile_plugin_t* operator()() const
            {
            return (molfile_plugin_t*)plugin_;
            }

        //! Destructor
        ~GSDPlugin()
            {
            fini();
            lib.close();
            }

    private:
        vmdplugin_t* plugin_; //!< VMD plugin
        DynamicLibrary lib; //!< Dynamic library for GSD plugin
        initfunc init_; //!< Function pointer to initialize the plugin
        finifunc fini_; //!< Function pointer to finalize the plugin

        //! Initialize the plugin
        bool init()
            {
            if(!lib.open()) return false;

            if (!init_)
                {
                // initialize the plugin library
                void *ifunc = lib.load("vmdplugin_init");
                if (ifunc)
                    {
                    init_ = (initfunc)ifunc;
                    }
                else
                    {
                    return false;
                    }
                }

            int initval = init_();
            return (initval == MOLFILE_SUCCESS);
            }

        //! Callback to set the plugin
        /*!
         * Must be static by design of the vmd api
         */
        static int setPlugin(void *v, vmdplugin_t *plugin)
            {
            GSDPlugin *self = (GSDPlugin *)v;
            self->plugin_ = plugin;
            return 0;
            }

        //! Finalize the plugin
        void fini()
            {
            if (!lib.open()) return;

            if (!fini_)
                {
                void *ffunc = lib.load("vmdplugin_fini");
                if (ffunc)
                    {
                    fini_ = (finifunc)ffunc;
                    }
                else
                    {
                    return;
                    }
                }
            fini_();
            }
    };

//! Test the dynamic library loader for the testing harness
TEST_CASE("DynamicLibrary")
    {
    SECTION("gsdplugin")
        {
        DynamicLibrary lib("./gsdplugin.so");

        // open the plugin library
        REQUIRE(lib.open());

        // initialize the plugin library
        void *ifunc = lib.load("vmdplugin_init");
        REQUIRE(ifunc);

        void *rfunc = lib.load("vmdplugin_register");
        REQUIRE(rfunc);

        void *ffunc = lib.load("vmdplugin_fini");
        REQUIRE(ffunc);

        lib.close();
        }

    SECTION("foobar")
        {
        DynamicLibrary lib("./foobar.so");

        // open the plugin library
        REQUIRE(!lib.open());

        // initialize the plugin library
        void *ifunc = lib.load("vmdplugin_init");
        REQUIRE(!ifunc);

        void *rfunc = lib.load("vmdplugin_register");
        REQUIRE(!rfunc);

        void *ffunc = lib.load("vmdplugin_fini");
        REQUIRE(!ffunc);

        lib.close();
        }
    }

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
    CHECK(plugin->minorv == 0);
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
        molfile_atom_t *atoms=(molfile_atom_t *)malloc(natoms*sizeof(molfile_atom_t));
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

        free(atoms);
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
    }
