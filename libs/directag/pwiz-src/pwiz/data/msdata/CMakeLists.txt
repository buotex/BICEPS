include_directories(${Boost_INCLUDE_DIRS})
include_directories(.)
include_directories("${PROJECT_SOURCE_DIR}/libs/directag/pwiz-src/pwiz")
SET(Directag_MSDATA_SRCS
    BinaryDataEncoder.cpp
    BinaryDataEncoderTest.cpp
    CVParam.cpp
    CVParamTest.cpp
    CVTranslator.cpp
    CVTranslatorTest.cpp
    cv.cpp
    ChromatogramList_mzML.cpp
    ChromatogramList_mzML_Test.cpp
    DefaultReaderList.cpp
    Diff.cpp
    DiffTest.cpp
    IO.cpp
    IOTest.cpp
    LegacyAdapter.cpp
    LegacyAdapterTest.cpp
    MSData.cpp
    MSDataFile.cpp
    RAMPAdapter.cpp
    RAMPAdapterTest.cpp
    Reader.cpp
    ReaderTest.cpp
    References.cpp
    ReferencesTest.cpp
    Serializer_mzML.cpp
    Serializer_mzML_Test.cpp
    Serializer_mzXML.cpp
    Serializer_mzXML_Test.cpp
    SpectrumIterator.cpp
    SpectrumIteratorTest.cpp
    SpectrumListFilter.cpp
    SpectrumListFilterTest.cpp
    SpectrumList_mzML.cpp
    SpectrumList_mzXML.cpp
    SpectrumList_mzXML_Test.cpp
)


#add_library(msdata MSData.cpp)
#add_library(reader Reader.cpp)
#add_library(spectrumiterator SpectrumIterator.cpp)


#file(GLOB folder_source *.cpp)
#foreach(source_file ${folder_source})
#  get_filename_component(filename ${source_file} NAME_WE)
#  add_library(${filename} ${source_file})
#  endforeach(source_file ${folder_source})
add_definitions(-DBOOST_IOSTREAMS_NO_LIB)
add_library(Directag_MSDATA ${Directag_MSDATA_SRCS})
