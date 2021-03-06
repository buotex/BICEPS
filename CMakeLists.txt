cmake_minimum_required(VERSION 2.4.7)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0015 NEW)
project(BICEPS)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(CMAKE_IS_EXPERIMENTAL YES_I_KNOW)
find_package(Boost COMPONENTS filesystem system regex program_options iostreams serialization thread REQUIRED)
#############################################################################
# build libraries
#############################################################################


OPTION(HAVE_PEPNOVO "ON, if you want to install and use Pepnovo"
	ON)
OPTION(HAVE_DIRECTAG "ON, if you want to install and use Directag"
	ON)


INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD)
IF(HAVE_UNISTD)
	ADD_DEFINITIONS(-DHAVE_UNISTD_H=${HAVE_UNISTD})
ENDIF(HAVE_UNISTD)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


#CHECK_INCLUDE_FILES(tcmalloc.h HAVE_TCMALLOC)
#IF(HAVE_TCMALLOC)
#    SET(TCMALLOC_LIB tcmalloc)
#ENDIF(HAVE_TCMALLOC)


INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(memmove HAVE_MEMMOVE)
IF(HAVE_MEMMOVE)
	ADD_DEFINITIONS(-DHAVE_MEMMOVE=${HAVE_MEMMOVE})
ENDIF(HAVE_MEMMOVE)


if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING
		"Choose the type of build, options are: Debug Release
		RelWithDebInfo MinSizeRel."
		FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# more warnings against everything, setting compiler-dependent flags
if(NOT MSVC)
	SET(CMAKE_CXX0X_FLAGS "-std=c++0x" )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare -Wno-sign-conversion " )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX0X_FLAGS}" )
	SET(BICEPS_CONFIG_PATH "getenv(\"HOME\")")
	SET(BICEPS_CONFINSTALL_PATH "$ENV{HOME}/.biceps")
ELSE(NOT MSVC)

	SET(BICEPS_CONFINSTALL_PATH "bin/.biceps")


endif(NOT MSVC)



CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/include/bicepsdefinitions.h.cmake
	${BICEPS_BINARY_DIR}/include/bicepsdefinitions.h
	)

add_subdirectory(libs)
add_subdirectory(src)


#############################################################################
# documentation
#############################################################################

# doxygen support
CONFIGURE_FILE(${BICEPS_SOURCE_DIR}/doc/Doxyfile.cmake
	${BICEPS_BINARY_DIR}/doc/Doxyfile
	@ONLY IMMEDIATE
	)
FIND_PACKAGE(Doxygen)

# target for doxygen
IF(DOXYGEN_FOUND)
	ADD_CUSTOM_TARGET(doc COMMAND ${DOXYGEN} ${BICEPS_BINARY_DIR}/doc/Doxyfile VERBATIM)
ENDIF(DOXYGEN_FOUND)



#############################################################################
# installation
#############################################################################



# documentation
INSTALL(DIRECTORY ${BICEPS_BINARY_DIR}/doc
	DESTINATION share/biceps
	COMPONENT documentation
	PATTERN "Doxyfile" EXCLUDE)
INSTALL(DIRECTORY ${BICEPS_SOURCE_DIR}/config/Models
	DESTINATION ${BICEPS_CONFINSTALL_PATH}
	COMPONENT config)
INSTALL(FILES ${BICEPS_SOURCE_DIR}/config/in_AAmodifications.param
	DESTINATION ${BICEPS_CONFINSTALL_PATH}
	COMPONENT config)


#############################################################################
# syncing docs
#############################################################################
ADD_CUSTOM_TARGET(syncdocs COMMAND ${BICEPS_SOURCE_DIR}/scripts/sync-docs-gh-pages.sh ${BICEPS_SOURCE_DIR} ${BICEPS_BINARY_DIR})
ADD_DEPENDENCIES(syncdocs doc)

#############################################################################
# packaging
#############################################################################


set(CPACK_PACKAGE_NAME "biceps")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A Novel Error Tolerant Search Strategy for Cross-Species Proteomics")


# CPack version strings

SET(CPACK_NSIS_MODIFY_PATH ON)
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Biceps")
SET(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})

SET(CPACK_RESOURCE_FILE_README "${BICEPS_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${BICEPS_SOURCE_DIR}/LICENSE.txt")

set(CPACK_COMPONENTS_ALL libraries headers applications config documentation)
set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "Application")
set(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "C++ Interface Documentation")

set(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION
	"The BICEPS executable.")
set(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION
	"The BICEPS documentation (Doxygen).")

INCLUDE(CPack)

