SET(
	DIRECTAG_SRC_DIRS
	${DIRECTAG_SOURCE_DIR}/freicore 
	${DIRECTAG_SOURCE_DIR}/pwiz-src
    ${DIRECTAG_SOURCE_DIR}/pwiz-src/libraries
    ${DIRECTAG_SOURCE_DIR}/pwiz-src/pwiz
	${DIRECTAG_SOURCE_DIR}/ext/src/expat-2.0.1/lib
	)
#include_directories(${DIRECTAG_SRC_DIR}/freicore ${DIRECTAG_SRC_DIR}/pwiz-src ${DIRECTAG_SRC_DIR}/pwiz-src/libraries)
include_directories(${DIRECTAG_SRC_DIRS})
link_directories(${DIRECTAG_BINARY_DIR}/freicore ${DIRECTAG_BINARY_DIR}/pwiz-src ${DIRECTAG_BINARY_DIR}/pwiz-src/libraries ${DIRECTAG_BINARY_DIR}/pwiz-src/libraries)
#add_executable(directag directag.cpp directagSpectrum.cpp)
add_library(directag directag.cpp directagSpectrum.cpp)

