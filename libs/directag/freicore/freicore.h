#ifndef _FREICORE_H
#define _FREICORE_H

#include "Profiler.h"
#include "ResidueMap.h"
#include "lnFactorialTable.h"
#include "shared_types.h"
#include "shared_defs.h"
#include "shared_funcs.h"
#include "SearchSpectrum.h"
#include "proteinStore.h"
#include "BaseSpectrum.h"

//#define BOOST_LIB_DIAGNOSTIC

#ifdef USE_MPI
	#undef SEEK_SET
	#undef SEEK_CUR
	#undef SEEK_END
	#include "mpi.h"
#endif

#endif
