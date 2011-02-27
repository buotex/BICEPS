PROJECT(EXPAT) # source files for zlib 
SET(EXPAT_SRCS xmlparse.c xmlrole.c xmltok.c xmltok_impl.c xmltok_ns.c ) 
ADD_LIBRARY(EXPAT ${EXPAT_SRCS})