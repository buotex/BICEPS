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
#target_link_libraries(directag freicore boost_filesystem boost_regex boost_program_options boost_iostreams boost_serialization boost_system Reader MSDataFile MSData DefaultReaderList vendor_readers IntegerSet CVParam cv Serializer_mzXML Serializer_mzML minixml LegacyAdapter pthread IO SpectrumList_mzML SpectrumList_mzXML ChromatogramList_mzML References BinaryDataEncoder zlib Base64 CVTranslator SHA1Calculator SHA1 spectrum_processing Diff z)

#target_link_libraries(directag freicore Reader MSDataFile MSData DefaultReaderList vendor_readers IntegerSet CVParam cv Serializer_mzXML Serializer_mzML minixml LegacyAdapter pthread IO SpectrumList_mzML SpectrumList_mzXML ChromatogramList_mzML References BinaryDataEncoder zlib Base64 CVTranslator SHA1Calculator SHA1 spectrum_processing Diff z)
#target_link_libraries(directag freicore Directag_MSDATA PWIZ_MISC_UTILITY)
