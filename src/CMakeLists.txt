
SET(
    DIRECTAG_SRC_DIRS
    ${BICEPS_SOURCE_DIR}/libs/directag/freicore 
    ${BICEPS_SOURCE_DIR}/libs/directag/pwiz-src 
    ${BICEPS_SOURCE_DIR}/libs/directag/pwiz-src/libraries
    ${BICEPS_SOURCE_DIR}/libs/directag/pwiz-src/pwiz
    ${BICEPS_SOURCE_DIR}/libs/directag/ext/src/expat-2.0.1/lib
   )
include_directories("${BICEPS_SOURCE_DIR}/libs/pepsplice" "${BICEPS_SOURCE_DIR}/libs/pepnovo" "${BICEPS_SOURCE_DIR}/libs/directag/directag" "${BICEPS_SOURCE_DIR}/libs/boost" ${DIRECTAG_SRC_DIRS} "${BICEPS_SOURCE_DIR}/libs/gmm-bic")
include_directories("${BICEPS_SOURCE_DIR}/include" "${BICEPS_BINARY_DIR}/include")
include_directories(${Boost_INCLUDE_DIRS})
  link_directories (${BICEPS_SOURCE_DIR}/src "${BICEPS_SOURCE_DIR}/libs/pepsplice" "${BICEPS_SOURCE_DIR}/libs/pepnovo" "${BICEPS_SOURCE_DIR}/libs/directag/directag" "${BICEPS_SOURCE_DIR}/libs" ${DIRECTAG_SRC_DIRS} "${BICEPS_SOURCE_DIR}/libs/gmm-bic")
#add_library(ets main.cpp fasta.cpp fasta.h)

if(MSVC)
  SET(
      WIN32_ONLY_LIBS
      thermo
     )
ELSE(MSVC)
  SET(
      WIN32_ONLY_LIBS
     )
  SET(
      PTHREAD
      pthread
     )

ENDIF(MSVC)

  IF(HAVE_PEPNOVO)
ADD_DEFINITIONS(-DHAVE_PEPNOVO)
  SET(PEPNOVO_LIBS
      pepnovo
     )
  ELSE(HAVE_PEPNOVO)
SET(PEPNOVO_LIBS)

ENDIF(HAVE_PEPNOVO)


  IF(HAVE_DIRECTAG)
ADD_DEFINITIONS(-DHAVE_DIRECTAG)
  SET(DIRECTAG_LIBS
      directag freicore Directag_MSDATA PWIZ_MISC_UTILITY vendor_readers minixml spectrum_processing ${ZLIB_SELECTOR} ${PTHREAD} ${WIN32_ONLY_LIBS}
      ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_REGEX_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
     )
  ELSE(HAVE_DIRECTAG)
  SET(DIRECTAG_LIBS)
ENDIF(HAVE_DIRECTAG)


  SET(GMMBIC_LIBS
      gmm_bic)



add_executable(biceps main.cpp biceps.cpp biceps.h fasta.cpp fasta.h commandparser.h commandparser.cpp)


target_link_libraries(biceps pepsplice ${Boost_PROGRAM_OPTIONS_LIBRARY} ${PEPNOVO_LIBS} ${DIRECTAG_LIBS} ${GMMBIC_LIBS}) 


  INSTALL(TARGETS biceps
      RUNTIME DESTINATION bin
      COMPONENT applications
      )

