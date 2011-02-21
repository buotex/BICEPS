include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}/libs/iostreams/src")
include_directories("${PROJECT_SOURCE_DIR}/libs/zlib")
ADD_LIBRARY(boost-zlib zlibdummy.cpp)
