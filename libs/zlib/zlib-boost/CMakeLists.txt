include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}/libs/iostreams/src")
include_directories("${BICEPS_SOURCE_DIR}/libs/zlib")
ADD_LIBRARY(zlib-boost zlibdummy.cpp)
