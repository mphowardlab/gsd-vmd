// Copyright (c) 2017, Michael P. Howard
// This file is part of the gsd-vmd project, released under the Modified BSD License.

#include "GSDPlugin.h"

GSDPlugin::GSDPlugin()
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

bool GSDPlugin::init()
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
int GSDPlugin::setPlugin(void *v, vmdplugin_t *plugin)
    {
    if (!v || !plugin) return 1;

    GSDPlugin *self = (GSDPlugin *)v;
    self->plugin_ = plugin;
    return 0;
    }

//! Finalize the plugin
void GSDPlugin::fini()
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
