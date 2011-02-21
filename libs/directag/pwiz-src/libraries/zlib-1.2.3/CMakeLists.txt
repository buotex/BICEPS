PROJECT(ZLIB) # source files for zlib 

ADD_SUBDIRECTORY(zlib-boost)
SET(ZLIB_SRCS adler32.c   gzio.c       inftrees.c  uncompr.c compress.c minigzip.c zutil.c crc32.c   deflate.c   inffast.c   inflate.c   trees.c ) 
ADD_LIBRARY(zlib ${ZLIB_SRCS})
