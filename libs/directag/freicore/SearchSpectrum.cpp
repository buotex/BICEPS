#include "stdafx.h"
#include "SearchSpectrum.h"
#include "searchResult.h"
#include "Profiler.h"

using namespace freicore;

SearchScoreList BaseSearchResult::emptyScoreList;

namespace std
{
	//ostream& operator<< ( ostream& o, const SearchScoreInfo& rhs )
	//{
	//	return ( o << static_cast< pair< string, float > >( rhs ) );
	//}

	//ostream& operator<< ( ostream& o, const SearchScoreList& rhs )
	//{
	//	return ( o << static_cast< vector< pair< string, float > > >( rhs ) );
	//}

	ostream& operator<< ( ostream& o, const BaseSearchResult& rhs )
	{
		return ( o << "(" << rhs.sequence << " " << rhs.getScoreList() << " " << rhs.key << " " << rhs.lociByIndex << " " << rhs.lociByName << ")" );
	}
}
