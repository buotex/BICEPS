#add_library(Base64 Base64.cpp)
#add_library(IntegerSet IntegerSet.cpp)
#add_library(Filesystem Filesystem.cpp)
#add_library(SHA1 SHA1.cpp)
#add_library(SHA1Calculator SHA1Calculator.cpp)
set(PWIZ_MISC_UTILITY_SRC
	Base64.cpp
	IntegerSet.cpp
	Filesystem.cpp
	SHA1.cpp
	SHA1Calculator.cpp
)
add_library(PWIZ_MISC_UTILITY ${PWIZ_MISC_UTILITY_SRC})