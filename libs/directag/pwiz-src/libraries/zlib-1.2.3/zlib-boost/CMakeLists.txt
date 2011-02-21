include_directories("${Boost_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIRS}/libs/iostreams/src")
include_directories("${PROJECT_SOURCE_DIR}/libs/directag/pwiz-src/libraries/zlib-1.2.3/zlib")
ADD_LIBRARY(zlib-boost zlibdummy.cpp)
