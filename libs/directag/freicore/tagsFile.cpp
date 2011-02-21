#include "stdafx.h"
#include "shared_defs.h"
#include "shared_funcs.h"
#include "tagsFile.h"

using namespace freicore;

namespace std
{
	ostream& operator<< ( ostream& o, const TagInfo& rhs )
	{
		o	<< "( "
			<< rhs.tag << ' '
			<< rhs.lowPeakMz << ' '
			<< rhs.nTerminusMass << ' '
			<< rhs.cTerminusMass << ' '
			<< rhs.valid << ' '
			<< rhs.totalScore;

		for( map< string, float >::const_iterator itr = rhs.scores.begin(); itr != rhs.scores.end(); ++itr )
			if( itr->first != "total" )
				o << ' ' << itr->second;
		return o << " )";
	}
}
