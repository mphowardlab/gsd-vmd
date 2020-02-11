// Copyright (c) 2017-2020, Michael P. Howard
// This file is part of the gsd-vmd project, released under the Modified BSD License.

#ifndef GSDVMD_TESTS_GSD_PLUGIN_H_
#define GSDVMD_TESTS_GSD_PLUGIN_H_

#include "DynamicLibrary.h"
#include "molfile_plugin.h"

extern "C"
    {
    //! VMD plugin init function
    typedef int (*initfunc)(void);
    //! VMD plugin registration function
    typedef int (*regfunc)(void *, vmdplugin_register_cb);
    //! VMD plugin fini function
    typedef int (*finifunc)(void);
    }

//! GSD plugin wrapper
class GSDPlugin
    {
    public:
        //! Constructor
        GSDPlugin();

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
        bool init();

        //! Callback to set the plugin
        /*!
         * Must be static by design of the vmd api
         */
        static int setPlugin(void *v, vmdplugin_t *plugin);

        //! Finalize the plugin
        void fini();
    };

#endif // GSDVMD_TESTS_GSD_PLUGIN_H_
