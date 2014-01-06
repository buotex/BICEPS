#ifndef _BASESPECTRUM_H
#define _BASESPECTRUM_H

#include "stdafx.h"
#include "base64.h"
#include "shared_defs.h"

using namespace freicore;

namespace freicore
{
	struct BaseSpectrum
	{
								BaseSpectrum();
								BaseSpectrum( const BaseSpectrum& old );

		virtual					~BaseSpectrum();

		float					CalculateComplementMz( float mz, int z );
		virtual size_t			size();

		SpectrumId				id;						    // source.index.charge
        string                  stringID;                   // arbitrary string identifier
        string                  nativeID;                   // exact native identifier
		string					fileName;					// name of file the spectrum was read from

		int						peakPreCount;				// number of peaks prior to preprocessing
		int						peakCount;					// number of peaks after preprocessing

		template< class Archive >
		void serialize( Archive& ar, const unsigned int version )
		{
			ar & id;
            ar & stringID;
            ar & nativeID;
			ar & fileName;
			ar & peakPreCount;
			ar & peakCount;
			ar & numSequenceComparisons;
			ar & numFragmentChargeStates;
			ar & mzOfPrecursor;
			ar & mOfPrecursor;
			ar & mOfUnadjustedPrecursor;
			ar & retentionTime;
			ar & processingTime;
			ar & mzUpperBound;
			ar & mzLowerBound;
			ar & totalIonCurrent;
			ar & totalPeakSpace;
		}

		int						numSequenceComparisons;		// number of times this spectrum has been compared to a sequence
		int						numFragmentChargeStates;

		float					mzOfPrecursor;				// mass/charge ratio
		float					mOfPrecursor;				// estimated neutral mass
		float					mOfUnadjustedPrecursor;		// calculated (fixed) neutral mass
		float					retentionTime;				// time in minutes this spectrum was seen
		float					processingTime;				// time in seconds this spectrum has taken to process
		float					mzUpperBound;				// the last peak in the spectrum before preprocessing
		float					mzLowerBound;				// the first peak in the spectrum before preprocessing
		float					totalIonCurrent;			// sum of absolute intensities of all peaks before preprocessing
		float					totalPeakSpace;				// the space between the first and last peak before preprocessing

	};

	struct spectraSortByOriginalPeakCount
	{
		bool operator() ( const BaseSpectrum* a, const BaseSpectrum* b )
		{
			return a->peakPreCount < b->peakPreCount;
		}
	};

	struct spectraSortByFilteredPeakCount
	{
		bool operator() ( const BaseSpectrum* a, const BaseSpectrum* b )
		{
			return a->peakCount < b->peakCount;
		}
	};

	struct spectraSortByID
	{
		bool operator() ( const BaseSpectrum* a, const BaseSpectrum* b )
		{
			if( a->id.source == b->id.source )
				if( a->id.index == b->id.index )
					return a->id.charge < b->id.charge;
				else
					return a->id.index < b->id.index;
			else
				return a->id.source < b->id.source;
		}
	};

	template< class SpectrumType, class SpectraListType >
	class BaseSpectraList : public list< SpectrumType* >
	{
	public:
		//typedef SpectraListType								ListType;
		typedef list< SpectrumType* >							BaseList;
		typedef typename BaseList::const_iterator				ListConstIterator;
		typedef typename BaseList::iterator						ListIterator;
		typedef map< SpectrumId, ListIterator >					ListIndex;
		typedef typename ListIndex::iterator					ListIndexIterator;

		ListIndex	index;

		BaseSpectraList() : BaseList() {}

		BaseSpectraList( const SpectraListType& rhs )
		{
			*this = rhs;
		}

		~BaseSpectraList()
		{
			clear();
		}

		SpectraListType& operator= ( const SpectraListType& rhs )
		{
			clear( false ); // if the spectra should be deleted, it should be before this point
			for( ListConstIterator itr = rhs.begin(); itr != rhs.end(); ++itr )
				BaseSpectraList< SpectrumType, SpectraListType >::push_back( *itr );
			//BaseList::insert( BaseList::end(), rhs.begin(), rhs.end() );
			return reinterpret_cast< SpectraListType& >( *this );
		}

		void insert( const ListIterator& startItr, const ListIterator& finishItr, const ListIterator& atItr )
		{
			//BaseList::insert( atItr, startItr, finishItr );
			for( ListIterator itr = startItr; itr != finishItr; ++itr )
				push_back( *itr );
		}

		void push_back( SpectrumType* s )
		{
			ListIndexIterator itr = index.find( s->id );
			if( itr != index.end() )
				cerr << "Warning: id \"" << s->id << "\" is already in the spectrum list" << endl;
			else
			{
				BaseList::push_back(s);
				index[ s->id ] = BaseList::end();
				-- index[ s->id ];
			}
		}

		void setId( const SpectrumId& oldId, const SpectrumId& newId )
		{
			if( oldId == newId )
				return;

			ListIndexIterator itr = index.find( oldId );
			if( itr != index.end() )
			{
				(*itr->second)->id = newId;
				index[ newId ] = itr->second;
				index.erase( itr );
			}
		}

		void erase( const ListIndexIterator& itr, bool deleteSpectrum = true )
		{
			erase( itr->second, deleteSpectrum );
		}

		void erase( const ListIterator& itr, bool deleteSpectrum = true )
		{
			ListIndexIterator indexItr = index.find( (*itr)->id );
			if( indexItr != index.end() )
				index.erase( indexItr );

			if( deleteSpectrum )
				delete *itr;

			BaseList::erase( itr );
		}

		void erase( const ListIterator& start, const ListIterator& end, bool deleteSpectrum = true )
		{
			for( ListIterator sItr = BaseList::begin(); sItr != BaseList::end(); ++sItr )
				erase( sItr++, deleteSpectrum );
		}

		void clear( bool deleteSpectra = true )
		{
			if( deleteSpectra )
				for( ListIterator sItr = BaseList::begin(); sItr != BaseList::end(); ++sItr )
					delete *sItr;

			index.clear();
			BaseList::clear();
		}

		void random_shuffle()
		{
			vector<SpectrumType*> v( this->begin(), this->end() );
			this->clear(false);
			std::random_shuffle( v.begin(), v.end() );
			for( typename vector<SpectrumType*>::const_iterator itr = v.begin(); itr != v.end(); ++itr )
				this->push_back( *itr );
		}

		void filterByChargeState(	int chargeState,
									SpectraListType* passingSpectra = NULL,
									SpectraListType* failingSpectra = NULL )
		{
			vector< int > chargeStates( 1, chargeState );
			filterByChargeState( chargeStates, passingSpectra, failingSpectra );
		}

		void filterByChargeState(	const vector< int >& chargeStates,
									SpectraListType* passingSpectra = NULL,
									SpectraListType* failingSpectra = NULL )
		{
			set<int> passByChargeState;
			for( size_t i=0; i < chargeStates.size(); ++i )
				passByChargeState.insert( chargeStates[i] );

			for( ListConstIterator sItr = BaseList::begin(); sItr != BaseList::end(); ++sItr )
			{
				SpectrumType* s = *sItr;

				if( passByChargeState.find( s->id.charge ) != passByChargeState.end() )
				{
					if( passingSpectra )
						passingSpectra->push_back( s );
				} else if( failingSpectra )
					failingSpectra->push_back( s );
			}
		}

		vector< size_t > getOriginalPeakCountStatistics()
		{
			vector< size_t > originalPeakCounts;
			for( ListIterator sItr = BaseList::begin(); sItr != BaseList::end(); ++sItr )
				originalPeakCounts.push_back( (size_t) (*sItr)->peakPreCount );

			std::sort( originalPeakCounts.begin(), originalPeakCounts.end() );

			vector< size_t > stats( 6 );
			stats[0] = originalPeakCounts.front();	// min
			stats[1] = originalPeakCounts.back();	// max
			stats[2] = originalPeakCounts[ originalPeakCounts.size() / 4 ];		// 1st quartile (25th percentile)
			stats[3] = originalPeakCounts[ originalPeakCounts.size() / 2 ];		// 2nd quartile (median)
			stats[4] = originalPeakCounts[ 3 * originalPeakCounts.size() / 4 ];	// 3rd quartile (75th percentile)
			stats[5] = accumulate( originalPeakCounts.begin(), originalPeakCounts.end(), 0 ) / originalPeakCounts.size(); // mean

			return stats;
		}

		vector< size_t > getFilteredPeakCountStatistics()
		{
			vector< size_t > filteredPeakCounts;
			for( ListIterator sItr = BaseList::begin(); sItr != BaseList::end(); ++sItr )
				filteredPeakCounts.push_back( (size_t) (*sItr)->peakCount );

			std::sort( filteredPeakCounts.begin(), filteredPeakCounts.end() );

			vector< size_t > stats( 6 );
			stats[0] = filteredPeakCounts.front();	// min
			stats[1] = filteredPeakCounts.back();	// max
			stats[2] = filteredPeakCounts[ filteredPeakCounts.size() / 4 ];		// 1st quartile (25th percentile)
			stats[3] = filteredPeakCounts[ filteredPeakCounts.size() / 2 ];		// 2nd quartile (median)
			stats[4] = filteredPeakCounts[ 3 * filteredPeakCounts.size() / 4 ];	// 3rd quartile (75th percentile)
			stats[5] = accumulate( filteredPeakCounts.begin(), filteredPeakCounts.end(), 0 ) / filteredPeakCounts.size(); // mean

			return stats;
		}
	};
}

// eliminate serialization overhead at the cost of never being able to increase the version.
BOOST_CLASS_IMPLEMENTATION( freicore::BaseSpectrum, boost::serialization::object_serializable )

// eliminate object tracking at the risk of a programming error creating duplicate objects.
BOOST_CLASS_TRACKING( freicore::BaseSpectrum, boost::serialization::track_never )

#endif
