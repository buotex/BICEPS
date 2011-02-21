#ifndef _SEARCHRESULT_H
#define _SEARCHRESULT_H

#include "shared_defs.h"
#include "shared_funcs.h"

namespace freicore
{
	typedef pair< string, float >		SearchScoreInfo;
	typedef vector< SearchScoreInfo >	SearchScoreList;
	struct BaseSearchResult
	{
		BaseSearchResult( const MvIntKey& k = MvIntKey(), const string& seq = "" )
			:	rank( 0 ), mass( 0.0f ), mod( 0.0f ), key( k ), sequence( seq )
		{}

		BaseSearchResult( const CandidateSequenceInfo& c )
			:	rank( 0 ), mass( c.mass ), sequence( c.sequence ),
				numTerminiCleavages( c.numTerminiCleavages ),
				numMissedCleavages( c.numMissedCleavages )
		{}

		virtual ~BaseSearchResult() {};

		size_t				rank;
		float				mass;
		float				mod;
		MvIntKey			key;
		string				sequence;
		ProteinLociByIndex	lociByIndex;
		ProteinLociByName	lociByName;
		int					numTerminiCleavages;
		int					numMissedCleavages;

		static SearchScoreList emptyScoreList;

		virtual float					getTotalScore() const = 0;
		virtual SearchScoreList			getScoreList() const = 0;

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & rank & mass & mod;
			ar & key;
			ar & sequence;
			ar & lociByIndex & lociByName;
			ar & numTerminiCleavages & numMissedCleavages;
		}
	};
}

namespace std
{
	//ostream& operator<< ( ostream& o, const SearchScoreInfo& rhs );
	//ostream& operator<< ( ostream& o, const SearchScoreList& rhs );
	ostream& operator<< ( ostream& o, const BaseSearchResult& rhs );
}

namespace freicore
{
	template< class SearchResultType >
	struct SearchResultSet : public topset< SearchResultType >
	{
		typedef SearchResultType SearchResult;
		typedef topset< SearchResultType > BaseSet;

		SearchResultSet( size_t MaxResults = 0 ) : topset< SearchResult >( MaxResults ) {}

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & boost::serialization::base_object< topset< SearchResult > >( *this );
		}

		void add( const SearchResult& newResult )
		{
			TemplateSetInsertPair( SearchResult ) rv = insert( newResult );
			if( rv.first != BaseSet::end() )
			{
				SearchResult& r = const_cast< SearchResult& >( *rv.first );
				r.lociByIndex.insert( newResult.lociByIndex.begin(), newResult.lociByIndex.end() );
				if( r.sequence.find( '+' ) != string::npos || newResult.sequence.find( '+' ) != string::npos )
				{
					SearchResult mergedResult( newResult );
					/*if( GetUnmodifiedSequence( r.sequence ) != GetUnmodifiedSequence( newResult.sequence ) )
						cout << r << endl << "should not combine with" << endl << newResult << endl << *this << endl;
					else*/
						mergedResult.sequence = MergeDeltaMasses( r.sequence, newResult.sequence );
					erase( rv.first );
					insert( mergedResult );
				}
			}
		}

		float getBestTotalScore() const
		{
			if( BaseSet::empty() )
				return 0;
			return const_cast< SearchResultType& >( *BaseSet::rbegin() ).getTotalScore();
		}

		string MergeDeltaMasses( const string& seq1, const string& seq2 )
		{
			string out;
			size_t maxLength = seq1.length() * 2;
			out.reserve( maxLength );
			size_t i1 = 0, i2 = 0;
			while( i1 < seq1.length() || i2 < seq2.length() )
			{
				if( i1 < seq1.length() && i2 < seq2.length() && seq1[i1] == seq2[i2] )
				{
					out.push_back( seq1[i1] );
					++i1; ++i2;
				} else
				{
					if( i1 < seq1.length() && seq1[i1] == '+' )
					{
						out.push_back( '+' );
						++i1;
					} else if( i2 < seq2.length() && seq2[i2] == '+' )
					{
						out.push_back( '+' );
						++i2;
					}
				}
			}
			return out;
		}

		void convertProteinIndexesToNames( const map< ProteinIndex, ProteinName >& indexToName )
		{
			for( typename BaseSet::iterator itr = BaseSet::begin(); itr != BaseSet::end(); ++itr )
				for( ProteinLociByIndex::iterator itr2 = itr->lociByIndex.begin(); itr2 != itr->lociByIndex.end(); ++itr2 )
				{
					ProteinName name;
					map< ProteinIndex, ProteinName >::const_iterator itr3 = indexToName.find( itr2->index );
					if( itr3 == indexToName.end() )
						name = lexical_cast<string>( itr2->index );
					else
						name = itr3->second;
					const_cast< SearchResult& >( *itr ).lociByName.insert( ProteinLocusByName( name, itr2->offset ) );
				}
		}

		void convertProteinNamesToIndexes( const map< ProteinName, ProteinIndex >& nameToIndex )
		{
			for( typename BaseSet::iterator itr = BaseSet::begin(); itr != BaseSet::end(); ++itr )
				for( ProteinLociByName::iterator itr2 = itr->lociByName.begin(); itr2 != itr->lociByName.end(); ++itr2 )
				{
					ProteinIndex index;
					map< ProteinName, ProteinIndex >::const_iterator itr3 = nameToIndex.find( itr2->name );
					if( itr3 == nameToIndex.end() )
						index = 0;
					else
						index = itr3->second;
					const_cast< SearchResult& >( *itr ).lociByIndex.insert( ProteinLocusByIndex( index, itr2->offset ) );
				}
		}

		void calculateRanks()
		{
			if( BaseSet::empty() )
				return;

			size_t rank = 1;
			float lastScore = BaseSet::rbegin()->getTotalScore();
			for( typename BaseSet::reverse_iterator itr = BaseSet::rbegin(); itr != BaseSet::rend(); ++itr )
			{
				const_cast< SearchResult& >( *itr ).rank = ( lastScore != itr->getTotalScore() ? ++rank : rank );
				lastScore = itr->getTotalScore();
			}
		}

		/*void calculateDeltaCn()
		{
			if( empty() )
				return;

			size_t numScores = scoreIndexMap.size();
			if( !scoreIndexMap.count("deltacn") )
				scoreIndexMap["deltacn"] = numScores;
			size_t deltaCnScoreIndex = scoreIndexMap["deltacn"];

			for( reverse_iterator itr = rbegin(); itr != rend(); ++itr )
			{
				SearchResult& r = const_cast< SearchResultT& >( *itr );
				if( r.scores.size() <= deltaCnScoreIndex )
					r.scores.push_back(0);
				r.scores[deltaCnScoreIndex] = ( r.scores[0] != 0 ? ( rbegin()->scores[0] - r.scores[0] ) / rbegin()->scores[0] : 0 );
			}
		}

		map< string, size_t > scoreIndexMap;*/
	};
}
BOOST_CLASS_IMPLEMENTATION( freicore::BaseSearchResult, boost::serialization::object_serializable )
//BOOST_CLASS_IMPLEMENTATION( freicore::SearchResultSet, boost::serialization::object_serializable )
BOOST_CLASS_TRACKING( freicore::BaseSearchResult, boost::serialization::track_never )
//BOOST_CLASS_TRACKING( freicore::SearchResultSet, boost::serialization::track_never )

namespace std
{
	template< class SearchResultType >
	ostream& operator<< ( ostream& o, const SearchResultSet< SearchResultType >& rhs )
	{
		return ( o << static_cast< set< SearchResultType > >( rhs ) );
	}
}

#endif
