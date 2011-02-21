#include "stdafx.h"

namespace freicore {
HostEndianType GetHostEndianType()
{
	int testInt = 127;
	char* testIntP = (char*) &testInt;

	if( testIntP[0] == 127 )
		return COMMON_LITTLE_ENDIAN;
	else if( testIntP[ sizeof(int)-1 ] == 127 )
		return COMMON_BIG_ENDIAN;
	else
		return COMMON_UNKNOWN_ENDIAN;
}

}
