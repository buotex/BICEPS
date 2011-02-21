#ifndef _DIRECTAG_H
#define _DIRECTAG_H
#include "stdafx.h"

#define DIRECTAG_MAJOR				1.
#define DIRECTAG_MINOR				2.
#define DIRECTAG_BUILD				2

#define DIRECTAG_VERSION				BOOST_PP_CAT( DIRECTAG_MAJOR, BOOST_PP_CAT( DIRECTAG_MINOR, DIRECTAG_BUILD ) )
#define DIRECTAG_VERSION_STRING			BOOST_PP_STRINGIZE( DIRECTAG_VERSION )
#define DIRECTAG_BUILD_DATE				"6/11/2008"

#define DIRECTAG_LICENSE				"Vanderbilt University Mass Spectrometry Research Center, D.Tabb/M.Chambers\n" \
							"Licensed under the Mozilla Public License.\n"

int directagFunc(int argc, vector<string> & argv, string & tags, vector<float> & lowPeakMzs, string & cachename);
#endif
