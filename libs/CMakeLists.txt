include_directories(${Boost_INCLUDE_DIRS}) 
add_subdirectory(pepsplice)
if(HAVE_PEPNOVO)
add_subdirectory(pepnovo)
endif(HAVE_PEPNOVO)
if(HAVE_DIRECTAG)
add_subdirectory(directag)
find_package(ZLIB)
if(NOT ZLIB_FOUND)
    add_subdirectory(zlib)
    set(ZLIB_SELECTOR ${ZLIB_SELECTOR} PARENT_SCOPE)
else(NOT ZLIB_FOUND)
    set(ZLIB_SELECTOR ${ZLIB_LIBRARIES} PARENT_SCOPE)
endif(NOT ZLIB_FOUND)

endif(HAVE_DIRECTAG)
add_subdirectory(gmm-bic)
