FILE(GLOB thermo_source *.cpp)
FILE(GLOB thermo_header *.h)
SOURCE_GROUP("Source Files" FILES ${thermo_source})
SOURCE_GROUP("Header Files" FILES ${thermo_header})

add_library(thermo ${thermo_source} ${thermo_header})