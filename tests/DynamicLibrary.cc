#include "DynamicLibrary.h"

#include <iostream>

/*!
 * \returns True if library is opened, false on failure
 */
bool DynamicLibrary::open()
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

/*!
 * \param name Symbol name
 * \returns Result of load on success, or nullptr on failure
 */
void* DynamicLibrary::load(const std::string& name) const
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

void DynamicLibrary::close()
    {
    if (handle_)
        {
        dlclose(handle_);
        handle_ = nullptr;
        }
    }