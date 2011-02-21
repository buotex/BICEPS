#ifndef EXPAT_XML_H
#define EXPAT_XML_H

#ifndef Expat_INCLUDED
	// MSVC magic to get Expat to link statically
	#if defined(WIN32) || defined(WIN64)
		#ifndef _DEBUG // Expat doesn't work in Debug mode... something about duplicate symbols
			#ifndef XML_STATIC
				#define XML_STATIC
			#endif
		#endif
	#endif

	#include "expat.h"
#endif

#endif
