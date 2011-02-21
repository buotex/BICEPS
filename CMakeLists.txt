cmake_minimum_required(VERSION 2.4.7)
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0015 NEW)
project(biceps)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(CMAKE_IS_EXPERIMENTAL YES_I_KNOW)
find_package(Boost COMPONENTS filesystem system regex program_options iostreams serialization thread REQUIRED)
if(!Boost_FOUND)
set(Boost_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/boost/boost)
set(Boost_LIBRARY_DIRS ${PROJECT_SOURCE_DIR}/lib)
endif(!Boost_FOUND)

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
find_package(tcmalloc)
link_directories(${tcmalloc_INCLUDE_DIRS})
set(TCMALLOC_LIB ${tcmalloc_LIBRARIES})


#CHECK_INCLUDE_FILES(tcmalloc.h HAVE_TCMALLOC)
#IF(HAVE_TCMALLOC)
#    SET(TCMALLOC_LIB tcmalloc)
#ENDIF(HAVE_TCMALLOC)


INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(memmove HAVE_MEMMOVE)
IF(HAVE_MEMMOVE)
	ADD_DEFINITIONS(-DHAVE_MEMMOVE=${HAVE_MEMMOVE})
ENDIF(HAVE_MEMMOVE)


include_directories(${Boost_INCLUDE_DIRS})
include_directories("${Boost_INCLUDE_DIRS}" "${Boost_LIBRARY_DIRS}" "${PROJECT_SOURCE_DIR}/src")
link_directories ("${Boost_LIBRARY_DIRS}" "${PROJECT_SOURCE_DIR}/pwiz-src/libraries/zlib-1.2.3")
add_subdirectory(libs)
add_subdirectory(src)
