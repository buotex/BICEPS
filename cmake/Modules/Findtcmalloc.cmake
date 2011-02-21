# -Try to look for google-Perftools

include(LibFindMacros)

libfind_pkg_check_modules(tcmalloc_PKGCONF tcmalloc)
find_library(tcmalloc_LIBRARY
    NAMES tcmalloc
    PATHS ${tcmalloc_PKGCONF_LIBRARY_DIRS}
)

set(tcmalloc_PROCESS_LIBS tcmalloc_LIBRARY)

libfind_process(tcmalloc)
