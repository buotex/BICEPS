include_directories("../pwiz-src/pwiz" "../ext" "../ext/src/expat-2.0.1/lib")
FILE(GLOB folder_source *.cpp)
FILE(GLOB folder_header *.h)
SOURCE_GROUP("Source Files" FILES ${folder_source})
SOURCE_GROUP("Header Files" FILES ${folder_header})

add_library(freicore ${folder_source} ${folder_header})